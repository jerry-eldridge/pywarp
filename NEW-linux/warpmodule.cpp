#include <Python.h>

#include <iostream>
using namespace std;

//#include <Windows.h>

#include <iomanip>
#include <list>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <exception>
//#pragma comment(lib, "python27.lib")

#define PI 3.1415926535897932384626433832795

#define _DEBUG_ false // set to true if needing to debug warp.
#define MSG(s) if (_DEBUG_) cout << "warp:" << s << endl;


//#define sprintf_s sprintf

int g_w = 900;
int g_h = 900;
bool N_8 = true; // toggle on/off 8-neighbor connectivity in the grid. N_8=false means N_4

int g_id = 0;

class Particle {
public:
	int id;
	double x[2];
	double v[2];
	double a[2];
	double F[2];
	double U;
	double mass;
	list<int> adj;
	Particle() { id = g_id++;  mass = 1; U = 0; }
	~Particle() { }
};

int g_M = 20;
class ScalarField {
public:
	int N;
	double *dat;
	ScalarField() { dat = 0; N = 0; }
	void alloc(int M) { N = M*M; dat = new double[N]; int i; for (i = 0; i < N; i++) dat[i] = 0; }
	double Get(int i, int j) { return dat[i + g_M*j]; }
	void Set(int i, int j, double val) { dat[i + g_M*j] = val; }
	void dealloc() { if (dat) delete[]dat; }
	~ScalarField() { dealloc(); }
	inline double Min() { int i; double mn=1e8; for (i = 0; i < N; i++) if (dat[i] <mn) mn = dat[i]; return mn; }
	inline double Max() { int i; double mx=-1e8; for (i = 0; i < N; i++) if (dat[i] >mx) mx = dat[i]; return mx; }
};

class WarpPoint {
public:
	double x;
	double y;
	WarpPoint() { x = y = 0; }
	WarpPoint(double xx, double yy) { x = xx; y = yy; }
	~WarpPoint() {}
};
class WarpScalar {
public:
	double val[3];
	WarpScalar() { val[0] = 0; val[1] = 0; val[2] = 0; }
	WarpScalar(double x1, double x2, double x3) { val[0] = x1; val[1] = x2; val[2] = x3; }
	~WarpScalar() {}
};

class WarpImage {
public:
	int w;
	int h;
	unsigned char *b;
	unsigned char *g;
	unsigned char *r;
	WarpImage() {  w = h = 0; b = g = r = 0; }
	WarpImage(int width, int height) {
		alloc(width, height);
	}
	~WarpImage() {
		dealloc();
	}
	void alloc(int width, int height)
	{
		w = width;
		h = height;
		b = new unsigned char[w*h]; g = new unsigned char[w*h]; r = new unsigned char[w*h];
		int i;
		for (i = 0; i < w*h; i++) {
			b[i] = g[i] = r[i] = (unsigned char)0;
		}
	}
	void dealloc() { 
		if (b) delete[]b; if (g) delete[]g; if (r) delete[]r; w = h = 0; }
	void Set(int i, int j, WarpScalar color) {
		int n = j + h*i;
		b[n] = static_cast<unsigned char>(color.val[0]);
		g[n] = static_cast<unsigned char>(color.val[1]);
		r[n] = static_cast<unsigned char>(color.val[2]);
	}
	void Get(int i, int j, WarpScalar &color) {
		int n = j + h*i;
		color.val[0] = double(int(b[n]));
		color.val[1] = double(int(g[n]));
		color.val[2] = double(int(r[n]));
	}
	void Fill(WarpScalar color) {
		int i, j;
		for (j = 0; j < h; j++) {
			for (i = 0; i < w; i++) {
				Set(i, j, color);
			}
		}
	}
};

double MapTo(double x1, double y1, double x2, double y2, double x);
void DrawRectangle(WarpImage *img, WarpPoint *p, WarpScalar color);
void InterpolateRect(WarpPoint *p, WarpPoint *q, WarpPoint *r, float t);
void Interp(WarpPoint p, WarpPoint q, WarpPoint r, float t);
void SetWarpPixel(WarpImage *src, WarpImage *dst, WarpPoint *p, WarpPoint *q, float u, float w);
WarpPoint Interp(WarpPoint p, WarpPoint q, float t);
void DrawTexture(WarpImage *src, WarpImage *dst, WarpPoint *p, WarpPoint *q);
float Distance(WarpPoint p, WarpPoint q);

void WarpClose();
void WarpInitialize(double dt, double *u, int N, int w, int h);
void WarpStep(double k_spr, bool grid, double dt);
void Show(char *window, double *img);
void Initialize(double *u, Particle *part, int N);
void Schwartzchild_U(double mass, double *u_dat, int N);
void Step(WarpImage *src2, WarpImage *warp, Particle *part, ScalarField *U,
	double k_spr, int g_w, int g_h, double dt, bool grid);
void InitializeParticles(Particle *part, ScalarField *U, int g_M, int g_w, int g_h);
void UpdateSpringForces(Particle *part, int N, double lo, double hi, double k_spr, int g_w, int g_h);
void Animate(Particle *part, int N, double dt2);
void Warp(Particle *part, int N, WarpImage *src, WarpImage *dst);

Particle *g_part = 0;
double t, dt;
WarpImage *g_warp = 0;
WarpImage *g_src2 = 0;
ScalarField g_U;
double g_dt = 0.0125;// 0.0025;// 0.005;
bool g_grid = false;
double g_k_spr = 40;
double g_t = 0;
char g_fn[4096] = "";
double g_speckles = 0.5;

double Schwartzschild_Curvature(double r, double theta, double M)
{
	double R;
	static double G = 6.673e-11;
	static double c = 299792458.0;
	R = -G*M*(-1.0*G*M + 0.5*pow(c,2.0) * r) / (pow(c,2.0) * pow(r,3.0) * (2 * G*M - pow(c,2.0) * r)) - (-0.5*G*M / (pow(c,2.0) * r) +
		1.0 + 0.25 / pow(tan(theta),2.0)) / pow(r,2.0) + (0.5*G*M*pow(sin(theta),2.0) * tan(theta) + 0.25*pow(c,2.0) *
		r*sin(theta)*cos(theta) - 0.5*pow(c,2.0) * r*pow(cos(theta),2.0) * tan(theta)) / (pow(c,2.0) * pow(r,3.0) * pow(sin(theta),2.0)
		* tan(theta)) + (-1.5*G*M + 0.5*pow(c,2.0) * r)*(2 * G*M - pow(c,2.0) * r) / (pow(c,2.0) * pow(r,3.0) * (-2.0*G*M + 1.0*pow(c,2.0) * r));
	return R;
}

double Min(double a, double b) { if (a < b) { return a; } else { return b; } }
double Max(double a, double b) { if (a > b) { return a; } else { return b; } }

extern "C" {

	static PyObject *
		warp_init(PyObject *self, PyObject *args)
	{
		MSG("warp_init");
		PyObject *Obj0 = PyTuple_GetItem(args, 0); // image data
		PyObject *Obj1 = PyTuple_GetItem(args, 1); // M
		PyObject *Obj2 = PyTuple_GetItem(args, 2); // width w
		PyObject *Obj3 = PyTuple_GetItem(args, 3); // height

		long M = PyLong_AsLong(Obj1);
		g_M = M;
		long w = PyLong_AsLong(Obj2);
		long h = PyLong_AsLong(Obj3);

		int i, j;
		g_src2 = new WarpImage;
		g_warp = new WarpImage;
		g_src2->alloc(w, h);
		Py_ssize_t N1 = PyList_Size(Obj0);
		assert(N1 == 3*w*h);
		double *b = new double[w*h];
		double *g = new double[w*h];
		double *r = new double[w*h];
		PyObject *pxl_b = 0;
		PyObject *pxl_g = 0;
		PyObject *pxl_r = 0;
		int idx = 0;
		for (j = 0; j<N1; j+=3) {
			pxl_b = PyList_GetItem(Obj0, j);
			pxl_g = PyList_GetItem(Obj0, j+1);
			pxl_r = PyList_GetItem(Obj0, j+2);
			b[idx] = PyFloat_AsDouble(pxl_b);
			g[idx] = PyFloat_AsDouble(pxl_g);
			r[idx] = PyFloat_AsDouble(pxl_r);
			idx++;
		}
		for (j = 0; j < h; j++) {
			for (i = 0; i < w; i++) {
				int n = i + w*j; // row-major read in
				WarpScalar color(b[n], g[n], r[n]);
				g_src2->Set(i, j, color);
			}
		}
		delete[]b;
		delete[]g;
		delete[]r;

		g_w = w;
		g_h = h;

		int N = g_M*g_M;
		double *u = new double[N];
		double M_sun = 1.989e30;
		double M_earth = 5.972e24;
		double M_milky_way_black_hole = 4.31e6*M_sun;
		Schwartzchild_U(M_milky_way_black_hole, u, N);
		WarpInitialize(g_dt, u, N, g_w, g_h);
		delete[]u;

		char buff[4096] = "";
/*
		sprintf_s(buff, 
			"Warp is Open: Read in: %s. Allocated Structures. M = %d, w = %d, h = %d", 
			g_fn, M, w, h);
*/
		return Py_BuildValue("s", buff);
	}
	static PyObject *
		warp_set_u_schwartzchild(PyObject *self, PyObject *args)
	{
		MSG("warp_set_u_schwartzchild");
		PyObject *Obj0 = PyTuple_GetItem(args, 0); // image data
		double mass = PyFloat_AsDouble(Obj0);
		int N = g_M*g_M;
		double *u = new double[N];
		Schwartzchild_U(mass, u, N);
		WarpInitialize(g_dt, u, N, g_w, g_h);
		delete[]u;
		return Py_BuildValue("i", N);
	}

	static PyObject *
		warp_close(PyObject *self, PyObject *args)
	{
		MSG("warp_close");
		WarpClose();

		char buff[4096] = "";
/*
		sprintf_s(buff, "Warp is now Closed. Deallocated Structures.");
*/
		return Py_BuildValue("s", buff);
	}
	static PyObject *
		warp_step(PyObject *self, PyObject *args)
	{
		MSG("warp_step");
		WarpStep(g_k_spr, g_grid, g_dt);

		char buff[4096] = "";
/*
		sprintf_s(buff, "Warp Step");
*/
		return Py_BuildValue("s", buff);
	}

	static PyObject *
		warp_set_u(PyObject *self, PyObject *args)
	{
		MSG("warp_set_u");
		PyObject *lst = PyTuple_GetItem(args, 0);
		Py_ssize_t N = PyList_Size(lst);
		double *u = new double[N];
		int j;
		for (j = 0; j<N; j++) {
			PyObject *item;
			item = PyList_GetItem(lst, j);
			u[j] = PyFloat_AsDouble(item);
		}

		WarpInitialize(g_dt, u, N, g_w, g_h);
		delete[]u;
		return Py_BuildValue("i", N);
	}

	static PyObject *
	warp_set_x(PyObject *self, PyObject *args)
	{
		MSG("warp_set_x");
		PyObject *lst = PyTuple_GetItem(args, 0);
		Py_ssize_t N = PyList_Size(lst);
		int j;
		for (j = 0; j<N; j++) {
			PyObject *item;
			item = PyList_GetItem(lst, j);
			PyObject *x0;
			PyObject *x1;
			x0 = PyList_GetItem(item, 0);
			x1 = PyList_GetItem(item, 1);
			g_part[j].x[0] = PyFloat_AsDouble(x0);
			g_part[j].x[1] = PyFloat_AsDouble(x1);

		}
		return Py_BuildValue("i", N);
	}

	static PyObject *
		warp_set_v(PyObject *self, PyObject *args)
	{
		MSG("warp_set_v");
		PyObject *lst = PyTuple_GetItem(args, 0);
		Py_ssize_t N = PyList_Size(lst);
		int j;
		for (j = 0; j<N; j++) {
			PyObject *item;
			item = PyList_GetItem(lst, j);
			PyObject *v0;
			PyObject *v1;
			v0 = PyList_GetItem(item, 0);
			v1 = PyList_GetItem(item, 1);
			g_part[j].v[0] = PyFloat_AsDouble(v0);
			g_part[j].v[1] = PyFloat_AsDouble(v1);

		}
		return Py_BuildValue("i", N);
	}

	static PyObject *
		warp_set_F(PyObject *self, PyObject *args)
	{
		MSG("warp_set_F");
		PyObject *lst = PyTuple_GetItem(args, 0);
		Py_ssize_t N = PyList_Size(lst);
		int j;
		for (j = 0; j<N; j++) {
			PyObject *item;
			item = PyList_GetItem(lst, j);
			PyObject *F0;
			PyObject *F1;
			F0 = PyList_GetItem(item, 0);
			F1 = PyList_GetItem(item, 1);
			g_part[j].F[0] = PyFloat_AsDouble(F0);
			g_part[j].F[1] = PyFloat_AsDouble(F1);

		}
		return Py_BuildValue("i", N);
	}

	static PyObject *
		warp_set_mass(PyObject *self, PyObject *args)
	{
		MSG("warp_set_mass");
		PyObject *lst = PyTuple_GetItem(args, 0);
		Py_ssize_t N = PyList_Size(lst);
		int j;
		for (j = 0; j<N; j++) {
			PyObject *item;
			item = PyList_GetItem(lst, j);
			g_part[j].mass = PyFloat_AsDouble(item);
		}
		return Py_BuildValue("i", N);
	}

	static PyObject *
		warp_set_k_spr(PyObject *self, PyObject *args)
	{
		MSG("warp_set_k_spr");
		PyObject *item = PyTuple_GetItem(args, 0);
		double k_spr = PyFloat_AsDouble(item);
		g_k_spr = k_spr;
		return Py_BuildValue("d", k_spr);
	}
	static PyObject *
		warp_set_dt(PyObject *self, PyObject *args)
	{
		MSG("warp_set_dt");
		PyObject *item = PyTuple_GetItem(args, 0);
		double dt = PyFloat_AsDouble(item);
		g_dt = dt;
		return Py_BuildValue("d", dt);
	}
	static PyObject *
		warp_get_t(PyObject *self, PyObject *args)
	{
		MSG("warp_get_t");
		return Py_BuildValue("d", g_t);
	}
	static PyObject *
		warp_get_width(PyObject *self, PyObject *args)
	{
		MSG("warp_get_width");
		return Py_BuildValue("i", g_w);
	}
	static PyObject *
		warp_get_height(PyObject *self, PyObject *args)
	{
		MSG("warp_get_height");
		return Py_BuildValue("i", g_h);
	}
	static PyObject *
		warp_get_u(PyObject *self, PyObject *args)
	{
		MSG("warp_get_u");
		int i, M, N;
		N = g_U.N;
		M = int(sqrt(N));
		PyObject *L = PyList_New(N);
		for (i = 0; i < N; i++) {
			PyList_SetItem(L, i, Py_BuildValue("d", g_U.dat[i]));
		}
		return L;
	}
	static PyObject *
		warp_get_x(PyObject *self, PyObject *args)
	{
		MSG("warp_get_x");
		int i, M, N;
		N = g_U.N;
		M = int(sqrt(N));
		PyObject *L = PyList_New(N);
		for (i = 0; i < N; i++) {
			PyObject *pt = PyList_New(2);
			PyList_SetItem(pt, 0, Py_BuildValue("d", g_part[i].x[0]));
			PyList_SetItem(pt, 1, Py_BuildValue("d", g_part[i].x[1]));
			PyList_SetItem(L, i, pt);
		}
		return L;
	}
	static PyObject *
		warp_get_v(PyObject *self, PyObject *args)
	{
		MSG("warp_get_v");
		int i, M, N;
		N = g_U.N;
		M = int(sqrt(N));
		PyObject *L = PyList_New(N);
		for (i = 0; i < N; i++) {
			PyObject *pt = PyList_New(2);
			PyList_SetItem(pt, 0, Py_BuildValue("d", g_part[i].v[0]));
			PyList_SetItem(pt, 1, Py_BuildValue("d", g_part[i].v[1]));
			PyList_SetItem(L, i, pt);
		}
		return L;
	}
	static PyObject *
		warp_get_F(PyObject *self, PyObject *args)
	{
		MSG("warp_get_F");
		int i, M, N;
		N = g_U.N;
		M = int(sqrt(N));
		PyObject *L = PyList_New(N);
		for (i = 0; i < N; i++) {
			PyObject *pt = PyList_New(2);
			PyList_SetItem(pt, 0, Py_BuildValue("d", g_part[i].F[0]));
			PyList_SetItem(pt, 1, Py_BuildValue("d", g_part[i].F[1]));
			PyList_SetItem(L, i, pt);
		}
		return L;
	}

	static PyObject *
		warp_get_image(PyObject *self, PyObject *args)
	{
		MSG("warp_get_image");
		int i, j, w, h;
		w = g_w;
		h = g_h;
		PyObject *L = PyList_New(w*h);
		for (j = 0; j < h; j++) {
			for (i = 0; i < w; i++) {
				WarpScalar val;
				g_warp->Get(i, j, val);
				double B, G, R;
				B = val.val[0];
				G = val.val[1];
				R = val.val[2];
				unsigned char b, g, r;
				b = (unsigned char)B;
				g = (unsigned char)G;
				r = (unsigned char)R;
				PyObject *pxl = PyList_New(3);
				PyList_SetItem(pxl, 0, Py_BuildValue("B", b));
				PyList_SetItem(pxl, 1, Py_BuildValue("B", g));
				PyList_SetItem(pxl, 2, Py_BuildValue("B", r));
				int n = i + w*j; // row major
				PyList_SetItem(L, n, pxl);
			}
		}
		return L;
	}

	static PyMethodDef WarpMethods[] = {
			{ "open", warp_init, METH_VARARGS,
			"open(fn, M, w, h). Open the warp allocating.\nfn is filename of image to warp.\nM is size of grid\nw x h image size. " },
			{ "close", warp_close, METH_VARARGS,
			"close(). Close the warp. Deallocates" },
			{ "step", warp_step, METH_VARARGS,
			"step(). Step Animate the warp. Warp and Display Grid." },
			{ "set_k_spr", warp_set_k_spr, METH_VARARGS,
			"set_k_spr(k_spr). Set Spring Constant k_spr." },
			{ "set_dt", warp_set_dt, METH_VARARGS,
			"set_dt(dt). Set Simulation Time Increment dt." },
			{ "get_t", warp_get_t, METH_VARARGS,
			"t = get_t(). Gets the simulation time t." },
			{ "get_width", warp_get_width, METH_VARARGS,
			"w = get_width(). Gets the image width." },
			{ "get_height", warp_get_height, METH_VARARGS,
			"h = get_height(). Gets the image height." },
			{ "set_u", warp_set_u, METH_VARARGS,
			"set_u(L). sets u with list L. Uses Initialization image for warping.\nGrid will warp based on scalar field u. nM = sqrt(len(L))" },
			{ "set_x", warp_set_x, METH_VARARGS,
			"set_x(L). sets x with list L = [[x0,x1],[x0,x1],...[x0,x1]] where len(L) = M*M." },
			{ "set_v", warp_set_v, METH_VARARGS,
			"set_v(L). sets v with list L = [[v0,v1],[v0,v1],...[v0,v1]] where len(L) = M*M." },
			{ "set_F", warp_set_F, METH_VARARGS,
			"set_F(L). sets F with list L = [[F0,F1],[F0,F1],...[F0,F1]] where len(L) = M*M." },
			{ "set_mass", warp_set_mass, METH_VARARGS,
			"set_mass(L). set mass with L = [m,m,m,...,m] where len(L) = M*M." },
			{ "set_u_schwartzchild", warp_set_u_schwartzchild, METH_VARARGS,
			"set_u_schwartzchild(mass_kg). sets u with the Schwartzchild metric's \nScalar Curvature given mass of spherical body in kg." },
			{ "get_u", warp_get_u, METH_VARARGS,
			"L = get_u(). gets the scalar field as a list L, size is M*M flattened." },
			{ "get_x", warp_get_x, METH_VARARGS,
			"L = get_x(). gets the positions L = [[x0,x1],[x0,x1],...[x0,x1]] where len(L) = M*M." },
			{ "get_v", warp_get_v, METH_VARARGS,
			"L = get_v). gets the velocities L = [[v0,v1],[v0,v1],...[v0,v1]] where len(L) = M*M." },
			{ "get_F", warp_get_F, METH_VARARGS,
			"L = get_F(). gets the forces L = [[F0,F1],[F0,F1],...[F0,F1]] where len(L) = M*M." },
			{ "image", warp_get_image, METH_VARARGS,
			"img = image(). gets warp image and returns a list img [[b,g,r],[b,g,r],...]." },
			{ NULL, NULL, 0, NULL }        /* Sentinel */
	};
        static struct PyModuleDef pywarp2 = {
             PyModuleDef_HEAD_INIT,
             "warp",
             "arg",
             -1,
             WarpMethods,
        };
	PyMODINIT_FUNC PyInit_warp(void)
	{
		return PyModule_Create(&pywarp2);
	}

}

void WarpClose()
{
	MSG("WarpClose()");
	delete[]g_part;

}

void WarpStep(double k_spr, bool grid, double dt)
{
	MSG("WarpStep");
	Step(g_src2, g_warp, g_part, &g_U, k_spr, g_w, g_h, dt, grid);
}

void WarpInitialize(double dt, double *u, int N, int w, int h)
{
	MSG("WarpInitialize");
	int M = int(sqrt(N));
	g_warp->dealloc();
	if (g_U.dat) g_U.dealloc();
	g_warp->alloc(w, h);

	g_U.alloc(M);
	g_part = new Particle[g_U.N];
	g_M = M;

	Initialize(u, g_part, g_U.N);
}

void Schwartzchild_U(double mass, double *u_dat, int N)
{
	MSG("Schwartzchild_U");
	int M = int(sqrt(N));
	int i, j;

	for (i = 0; i < M; i++) {
		double c = PI;
		double x = MapTo(0, -c, M - 1, c, i);
		for (j = 0; j < M; j++) {
			double y = MapTo(0, -c, M - 1, c, j);
			double r = sqrt(pow(x, 2.0) + pow(y, 2.0)) + 0.01; // make non-zero
			//double R = 0.5 / pow(r, 2.0);
			double theta = atan2(y, x);
			double R = abs(Schwartzschild_Curvature(r, theta, mass));
			double u = 1.0 / R;
			u = -log(u);
			u_dat[j + M*i] = u;
		}
	}

}

void Initialize(double *u, Particle *part, int N)
{
	MSG("Initialize");
	int i;
	for (i = 0; i < N; i++) g_U.dat[i] = u[i];

	InitializeParticles(part, &g_U, g_M, g_w, g_h);
}


// With x1 -> y1 and x2 -> y2, given x, return y using linear map
double MapTo(double x1, double y1, double x2, double y2, double x)
{
	//MSG("MapTo");
	double epsilon = 0.0001;
	//if (_DEBUG_) cout << x1 << "->" << y1 << "," << x2 << "->" << y2 << " : " << x << endl;
	try {
		if (abs(x2 - x1) > epsilon) {
			double m = 1.*(y2 - y1) / (x2 - x1);
			double y = m*(x - x1) + y1;
			return y;
		}
		else {
			double y = (y1 + y2) / 2.0;
			return y;
		}
	}
	catch (exception e)
	{
		cout << "Error: MapTo" << endl;
		system("pause");
	}
}


void DrawTexture(WarpImage *src, WarpImage *dst, WarpPoint *p, WarpPoint *q)
{
	//MSG("DrawTexture");
	float u, w, du, dw;
	// Optimize later du and dw to more appropriate size
	int dx, dy;
	float u1, w1;
	u1 = Distance(p[1], p[2]);
	u1 = Max(u1, Distance(p[0], p[3]));
	u1 = Max(u1, Distance(q[1], q[2]));
	u1 = Max(u1, Distance(q[0], q[3]));

	w1 = Distance(p[0], p[1]);
	w1 = Max(w1, Distance(p[3], p[2]));
	w1 = Max(w1, Distance(q[0], q[1]));
	w1 = Max(w1, Distance(q[3], q[2]));

	// Set g_speckles to 0.5
	du = g_speckles / (u1+1);
	dw = g_speckles / (w1+1);
	u = 0;
	while (u < 1) {
		w = 0;
		while (w < 1) {
			SetWarpPixel(src, dst, p, q, u, w);
			w += dw;
		}
		u += du;
	}
}

float Distance(WarpPoint p, WarpPoint q)
{
	//MSG("Distance");
	float dx, dy;
	dx = abs(p.x - q.x);
	dy = abs(p.y - q.y);
	//float dist = dx + dy;  // L1
	float dist = sqrt(dx*dx + dy*dy); // L2

	return dist;
}

void SetWarpPixel(WarpImage *src, WarpImage *dst, WarpPoint *p, WarpPoint *q, float u, float w)
{
	//MSG("SetWarpPixel");
	WarpPoint Q1 = Interp(q[0], q[3], u);
	WarpPoint Q2 = Interp(q[1], q[2], u);
	WarpPoint Q = Interp(Q2, Q1, w);

	WarpPoint P1 = Interp(p[0], p[3], u);
	WarpPoint P2 = Interp(p[1], p[2], u);
	WarpPoint P = Interp(P2, P1, w);

	// This longer method is supposed to be quicker than the commented out
	// Get2D and Set2D lines
	int i, j, ii, jj;
	i = P.x; j = P.y;
	ii = Q.x; jj = Q.y;

	if ((i >= src->w) || (i < 0)) return;
	if ((j >= src->h) || (j < 0)) return;

	if ((ii >= dst->w) || (ii < 0)) return;
	if ((jj >= dst->h) || (jj < 0)) return;

	WarpScalar color;
	src->Get(i, j, color);
	dst->Set(ii, jj, color);

}

void InterpolateRect(WarpPoint *p, WarpPoint *q, WarpPoint *r, float t)
{
	MSG("InterpolateRect");
	int i;
	for (i = 0; i<4; i++) {
		r[i] = Interp(p[i], q[i], t);
	}

}

WarpPoint Interp(WarpPoint p, WarpPoint q, float t)
{
	//MSG("Interp");
	WarpPoint r(0, 0);

	r.x = p.x*(1 - t) + q.x*t;
	r.y = p.y*(1 - t) + q.y*t;

	return r;
}

void InitializeParticles(Particle *part, ScalarField *U, int M, int g_w, int g_h) 
{
	MSG("InitializeParticles");
	int i;
	int N = M*M; 
	for (i = 0; i < N; i++) {
		part[i].x[0] = (i%M)*g_w / (M-1);
		part[i].x[1] = int(i / M)*g_h / (M-1);

		double x1 = part[i].x[0];
		double x2 = part[i].x[1];
		if (x1 == 0) x1 = 0.01;
		if (x2 == 0) x2 = 0.01;
		part[i].U = U->dat[i];

		part[i].v[0] = part[i].v[1] = 0;
		part[i].a[0] = part[i].a[0] = 0;
		if ((i + 1 < N) && (i%M != M - 1))  {
			part[i].adj.push_back(i + 1);
		}
		if ((i - 1 > 0) && (i%M != 0)) {
			part[i].adj.push_back(i - 1);
		}
		if (i + M < N)  {
			part[i].adj.push_back(i + M);
		}
		if (i - M >= 0) {
			part[i].adj.push_back(i - M);
		}

		if (N_8) {
			if ((i + M + 1 < N) && (i%M != M - 1))  {
				part[i].adj.push_back(i + M + 1);
			}
			if ((i + M - 1 < N) && (i%M != 0)) {
				part[i].adj.push_back(i + M - 1);
			}
			if ((i - (M - 1) >= 0) && (i%M != M - 1))  {
				part[i].adj.push_back(i - (M - 1));
			}
			if ((i - (M + 1) >= 0) && (i%M != 0)) {
				part[i].adj.push_back(i - (M + 1));
			}
		}

	}
}

void UpdateSpringForces(Particle *part, int N, double lo, double hi, double k_spr, int g_w, int g_h)
{
	MSG("UpdateSpringForces");
	int i;
	double Ax, Ay, Bx, By, Cx, Cy, norm, xx, len, Fx, Fy, u;
	int M;
	int n3;
	M = int(sqrt(N));
	//cout << "M = " << M << " ";
	// Calculate Spring Forces on Particles
	for (i = 0; i < N; i++) {
		list<int>::iterator qq;
		part[i].F[0] = part[i].F[1] = 0;
		for (qq = part[i].adj.begin(); qq != part[i].adj.end(); qq++) {
			n3 = *qq;
			//cout << "i = " << i << " ";
			//cout << "n3 = " << n3 << " ";
			Ax = part[i].x[0];
			Ay = part[i].x[1];
			Bx = part[n3].x[0];
			By = part[n3].x[1];
			Cx = Ax - Bx;
			Cy = Ay - By;
			norm = sqrt(Cx*Cx + Cy*Cy);
			//cout << "norm = " << norm << " ";
			double epsilon = 0.0001;
			u = (part[i].U + part[n3].U) / 2.0;
			//cout << "u = " << u << " ";
			if (abs(lo - hi) < epsilon) {
				len = 1;
			}
			else {
				len = MapTo(lo, g_w / M*01.45, hi, g_h / M*.55, u);
			}
			//cout << "len = " << len << " ";
			xx = norm - len;
			Fx = Fy = 0;
			Fx = -k_spr*xx*Cx / norm;
			Fy = -k_spr*xx*Cy / norm;
			//cout << "Fx, Fy = " << Fx << "," << Fy << endl;
			part[i].F[0] += Fx;
			part[i].F[1] += Fy;
		}
	}
	MSG("End of UpdateSpringForces")
}

void Animate(Particle *part, int N, double dt2)
{
	MSG("Animate");
	int k, i, j;
	for (i = 0; i < N; i++) {
		int c = part[i].adj.size();
		if (c < 8) continue;
		for (j = 0; j < 2; j++) {
			part[i].a[j] = part[i].F[j] / part[i].mass;
			part[i].v[j] = part[i].v[j] + part[i].a[j] * dt2;
			part[i].x[j] = part[i].x[j] + part[i].v[j] * dt2;
			part[i].v[j] = Min(10, Max(-10, part[i].v[j]));
			part[i].v[j] = 0.8*part[i].v[j];
		}
	}
}

void Warp(Particle *part, int N, WarpImage *src, WarpImage *dst)
{
	MSG("Warp");
	g_w = src->w;
	g_h = src->h;
	int i, j;
	int M = int(sqrt(N));
	WarpPoint p[4], q[4];
	for (i = 0; i < M; i++) {
		for (j = 0; j < M; j++) {
			int ii, k;
			ii = j + M*i; k = 0;
			p[k].x = (ii%M)*g_w / (M-1);
			p[k].y = int(ii / M)*g_h / (M-1);
			q[k].x = part[ii].x[0];
			q[k].y = part[ii].x[1];
			//cout << p[k].x << " " << p[k].y << " " << q[k].x << " " << q[k].y << endl;
			if ((i + 1) < M) {
				ii = j + M*(i + 1); k = 1;
				p[k].x = (ii%M)*g_w / (M-1);
				p[k].y = int(ii / M)*g_h / (M-1);
				q[k].x = part[ii].x[0];
				q[k].y = part[ii].x[1];
			}
			//cout << p[k].x << " " << p[k].y << " " << q[k].x << " " << q[k].y << endl;
			if ((j + 1) < M) {
				ii = j + 1 + M*i; k = 3;
				p[k].x = (ii%M)*g_w / (M - 1);
				p[k].y = int(ii / M)*g_h / (M - 1);
				q[k].x = part[ii].x[0];
				q[k].y = part[ii].x[1];
			}
			//cout << p[k].x << " " << p[k].y << " " << q[k].x << " " << q[k].y << endl;
			if ((j + 1) < M && (i + 1) < M) {
				ii = (j + 1) + M*(i + 1); k = 2;
				p[k].x = (ii%M)*g_w / (M - 1);
				p[k].y = int(ii / M)*g_h / (M - 1);
				q[k].x = part[ii].x[0];
				q[k].y = part[ii].x[1];
			}
			//cout << p[k].x << " " << p[k].y << " " << q[k].x << " " << q[k].y << endl;
			DrawTexture(src, dst, p, q);
		}
	}
}

void Step(WarpImage *src2, WarpImage *warp, Particle *part, ScalarField *U,
	double k_spr, int g_w, int g_h, double dt, bool grid)
{
	MSG("Step");
	g_t += dt;
	double uMin = U->Min();
	double uMax = U->Max();

	double R, G, B;
	R = 200; G = 200; B = 255;
	WarpScalar color(B, G, R);
	warp->Fill(color);
	UpdateSpringForces(part, U->N, uMin, uMax, k_spr, g_w, g_h);
	Animate(part, U->N, dt);
	Warp(part, U->N, src2, warp);
}
