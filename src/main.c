#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define EXPR_MAX_LENGTH 200

/* FUNCTION STRUCTURES */
typedef struct Polynomial {
	double x_coef;
	double x_exp;
} Poly;

typedef struct Exponential {
	double x_coef;
	double x_exp;
	double fn_coef;
	double fn_exp;
	double base;
} Exp;

typedef struct Logarithmic {
	double x_coef;
	double x_exp;
	double fn_coef;
	double fn_exp;
	double base;
} Log;

typedef struct Trigonometric {
	int trig_fn;
	double x_coef;
	double x_exp;
	double fn_coef;
	double fn_exp;
} Trig;

typedef struct Inverse_Trigonometric {
	int trig_fn;
	double x_coef;
	double x_exp;
	double fn_coef;
	double fn_exp;
} Inv_Trig;

typedef struct Function {
	Poly **poly_arr;
	int poly_count;
	Exp **exp_arr;
	int exp_count;
	Log **log_arr;
	int log_count;
	Trig **trig_arr;
	int trig_count;
	Inv_Trig **inv_trig_arr;
	int inv_trig_count;
} Function;

/* CONSTRUCTORS */
Poly *create_poly(double x_coef, double x_exp);
Exp *create_exp(double x_coef, double x_exp, double fn_coef, double fn_exp, double base);
Log *create_log(double x_coef, double x_exp, double fn_coef, double fn_exp, double base);
Trig *create_trig(int trig_fn, double x_coef, double x_exp, double fn_coef, double fn_exp);
Inv_Trig *create_inv_trig(int trig_fn, double x_coef, double x_exp, double fn_coef, double fn_exp);
Function *create_function(Poly **poly_arr, int poly_count, Exp **exp_arr, int exp_count, Log **log_arr,
	int log_count, Trig **trig_arr, int trig_count, Inv_Trig **inv_trig_arr, int inv_trig_count);
double **create_matrix(int n, int m);

/* INPUT FUNCTIONS */
void input_x(double *x_coef, double *x_exp);
void input_fn(double *fn_coef, double *fn_exp);
void input_poly(Poly **p);
void input_exp(Exp **e);
void input_log(Log **l);
void input_trig(Trig **t);
void input_inv_trig(Inv_Trig **it);
void input_function(Function **f);
void input_matrix(double ***matrix, int n, int m);

/* FUNCTION CALCULATORS */
double calculate_poly(Poly *p, double x);
double calculate_exp(Exp *e, double x);
double calculate_log(Log *l, double x);
double calculate_trig(Trig *t, double x);
double calculate_inv_trig(Inv_Trig *it, double x);
double calculate(Function *f, double x);

/* TO STRING FUNCTIONS */
char *to_string_poly(Poly *p);
char *to_string_exp(Exp *e);
char *to_string_log(Log *l);
char *to_string_trig(Trig *t);
char *to_string_inv_trig(Inv_Trig *it);
char *to_string_function(Function *f);

/* DESTRUCTORS */
void free_function(Function *f);
void free_matrix(double **matrix, int n);

/* PROCESSES */
double bisection(Function *f, double start, double end, double e, int opt, int max_iter);
double regula_falsi(Function *f, double start, double end, double e, int opt, int max_iter);
double newton_raphson(Function *f, double x0, double e, int max_iter);
double **inverse_matrix(double **matrix, int n);
double **gauss_elimination(double **a, double **c, int n);
double **gauss_seidel(double **a, double **c, double **starting_values, int n, double epsilon, int opt, int max_iter);
double backward_difference(Function *f, double x, double h);
double forward_difference(Function *f, double x, double h);
double centered_difference(Function *f, double x, double h);
double simpsons_rule_1_3(Function *f, double start, double end, int n);
double simpsons_rule_3_8(Function *f, double start, double end, int n);
double trapezoidal_rule(Function *f, double start, double end, int n);
Function *gregory_newton(int n, double x0, double h, double *vals);

/* HELPER FUNCTIONS */
double cot(double x);
double acot(double x);
void checkAllocation(void *ptr);
void print_matrix(double **matrix, int n, int m);
double determinant_matrix(double **matrix, int n);
double **minor_matrix(double **matrix, int n, int i, int j);
double **adjoint_matrix(double **matrix, int n);
double **augmented_matrix(double **a, double **c, int n);
int dominant_matrix_max_values(double **matrix, int n, int row, int col, double current, double *max, double *vals);
double **dominant_matrix(double **matrix, int n, int m, int opt);
Poly *multiply_poly(Poly *p1, Poly *p2);
Poly **multiply_poly_function(Poly **poly_arr1, int n1, Poly **poly_arr2, int n2);
Poly **merge_poly(Poly **poly_arr1, int n1, Poly **poly_arr2, int n2);
double fact(int n);

int main() {
	int choice;
	double x0, deltax, e, start, end, x, h;
	Function *f;
	int n, i, opt, max_iter;
	double **matrix, **matrix2, **matrix3, **matrix4;
	double *vals;
	char *str;

	choice = -1;
	while (choice != 0) {
		printf("\nQuit: 0\nBisection: 1\nRegula-Falsi: 2\nNewton Raphson: 3\nInverse Matrix: 4\nGauss Elimination: 5\nGauss-Seidel: 6\nNumerical Differentiation: 7\nSimpson's Rule: 8\nTrapezoidal Rule: 9\nGregory-Newton: 10\nChoice:\n");
		scanf(" %d", &choice);
		switch (choice) {
		case 1:
			input_function(&f);
			printf("start:\n");
			scanf(" %lf", &start);
			printf("end:\n");
			scanf(" %lf", &end);
			printf("epsilon:\n");
			scanf(" %lf", &e);
			opt = 0;
			while (opt < 1 || opt > 2) {
				printf("Stopping criterion:\nf(x) <= epsilon: 1\n(end - start) / 2^n <= epsilon: 2\nChoice:\n");
				scanf(" %d", &opt);
			}
			printf("Max iterations:\n");
			scanf("%d", &max_iter);
			printf("Result: %lf\n", bisection(f, start, end, e, opt, max_iter));
			free_function(f);
			break;
		case 2:
			input_function(&f);
			printf("start:\n");
			scanf(" %lf", &start);
			printf("end:\n");
			scanf(" %lf", &end);
			printf("epsilon:\n");
			scanf(" %lf", &e);
			opt = 0;
			while (opt < 1 || opt > 2) {
				printf("Stopping criterion:\nf(x) <= epsilon: 1\n(end - start) / 2^n <= epsilon: 2\nChoice:\n");
				scanf(" %d", &opt);
			}
			printf("Max iterations:\n");
			scanf("%d", &max_iter);
			printf("Result: %lf\n", regula_falsi(f, start, end, e, opt, max_iter));
			free_function(f);
			break;
		case 3:
			input_function(&f);
			printf("x0:\n");
			scanf(" %lf", &x0);
			printf("epsilon:\n");
			scanf(" %lf", &e);
			printf("Max iterations:\n");
			scanf("%d", &max_iter);
			printf("Result: %lf\n", newton_raphson(f, x0, e, max_iter));
			free_function(f);
			break;
		case 4:
			printf("N:\n");
			scanf(" %d", &n);
			input_matrix(&matrix, n, n);
			matrix2 = inverse_matrix(matrix, n);
			printf("Inverse Matrix: ");
			print_matrix(matrix2, n, n);
			free_matrix(matrix, n);
			free_matrix(matrix2, n);
			break;
		case 5:
			printf("\nA * X = C\n\n");
			printf("N:\n");
			scanf(" %d", &n);
			printf("A:\n");
			input_matrix(&matrix, n, n);
			printf("C:\n");
			input_matrix(&matrix2, n, 1);
			matrix3 = gauss_elimination(matrix, matrix2, n);
			printf("X: ");
			print_matrix(matrix3, n, 1);
			free_matrix(matrix, n);
			free_matrix(matrix2, n);
			free_matrix(matrix3, n);
			break;
		case 6:
			printf("\nA * X = C\n\n");
			printf("N:\n");
			scanf(" %d", &n);
			printf("A:\n");
			input_matrix(&matrix, n, n);
			printf("C:\n");
			input_matrix(&matrix2, n, 1);
			printf("Starting values:\n");
			input_matrix(&matrix3, n, 1);
			printf("epsilon:\n");
			scanf(" %lf", &e);
			printf("Find diagonally dominant matrix? 0/1:\n");
			scanf(" %d", &opt);
			printf("Max iterations:\n");
			scanf("%d", &max_iter);
			matrix4 = gauss_seidel(matrix, matrix2, matrix3, n, e, opt, max_iter);
			printf("X: ");
			print_matrix(matrix4, n, 1);
			free_matrix(matrix, n);
			free_matrix(matrix2, n);
			free_matrix(matrix3, n);
			free_matrix(matrix4, n);
			break;
		case 7:
			input_function(&f);
			opt = 0;
			while (opt < 1 || opt > 3) {
				printf("\nBackward Difference: 1\nForward Difference: 2\nCentered Difference: 3\nChoice:\n");
				scanf(" %d", &opt);
			}
			switch (opt) {
			case 1: printf("\n(f(x) - f(x - h)) / h\n"); break;
			case 2: printf("\n(f(x + h) - f(x)) / h\n"); break;
			case 3: printf("\n(f(x + h) - f(x - h)) / 2h\n"); break;
			}
			printf("x:\n");
			scanf(" %lf", &x);
			printf("h:\n");
			scanf(" %lf", &h);
			switch (opt) {
			case 1: printf("Result: %lf\n", backward_difference(f, x, h)); break;
			case 2: printf("Result: %lf\n", forward_difference(f, x, h)); break;
			case 3: printf("Result: %lf\n", centered_difference(f, x, h)); break;
			}
			free_function(f);
			break;
		case 8:
			input_function(&f);
			opt = 0;
			while (opt < 1 || opt > 2) {
				printf("\nSimpson's Rule 1/3: 1\nSimpson's Rule 3/8: 2\nChoice:\n");
				scanf(" %d", &opt);
			}
			printf("start:\n");
			scanf(" %lf", &start);
			printf("end:\n");
			scanf(" %lf", &end);
			printf("n:\n");
			scanf(" %d", &n);
			if (opt == 1) printf("Result: %lf\n", simpsons_rule_1_3(f, start, end, n));
			else printf("Result: %lf\n", simpsons_rule_3_8(f, start, end, n));
			free_function(f);
			break;
		case 9:
			input_function(&f);
			printf("start:\n");
			scanf(" %lf", &start);
			printf("end:\n");
			scanf(" %lf", &end);
			printf("n:\n");
			scanf(" %d", &n);
			printf("Result: %lf\n", trapezoidal_rule(f, start, end, n));
			free_function(f);
			break;
		case 10:
			printf("x0:\n");
			scanf(" %lf", &x0);
			printf("h:\n");
			scanf(" %lf", &h);
			printf("n:\n");
			scanf(" %d", &n);
			vals = (double *) malloc(n * sizeof(double));
			x = x0;
			for (i = 0; i < n; i++) {
				printf("f(%lf):\n", x);
				scanf(" %lf", &vals[i]);
				x += h;
			}
			f = gregory_newton(n, x0, h, vals);
			str = to_string_function(f);
			printf("\nResult: %s\n", str);
			printf("\nx:\n");
			scanf(" %lf", &x);
			printf("\nResult: %lf\n", calculate(f, x));
			free(vals);
			free_function(f);
			free(str);
			break;
		}
	}
	return 0;
}

/* CONSTRUCTORS */
Poly *create_poly(double x_coef, double x_exp) {
	Poly *p = (Poly *) malloc(sizeof(Poly));
	checkAllocation(p);
	p->x_coef = x_coef;
	p->x_exp = x_exp;
	return p;
}

Exp *create_exp(double x_coef, double x_exp, double fn_coef, double fn_exp, double base) {
	Exp *e = (Exp *) malloc(sizeof(Exp));
	checkAllocation(e);
	e->x_coef = x_coef;
	e->x_exp = x_exp;
	e->fn_coef = fn_coef;
	e->fn_exp = fn_exp;
	e->base = base;
	return e;
}

Log *create_log(double x_coef, double x_exp, double fn_coef, double fn_exp, double base) {
	Log *l = (Log *) malloc(sizeof(Log));
	checkAllocation(l);
	l->x_coef = x_coef;
	l->x_exp = x_exp;
	l->fn_coef = fn_coef;
	l->fn_exp = fn_exp;
	l->base = base;
	return l;
}

Trig *create_trig(int trig_fn, double x_coef, double x_exp, double fn_coef, double fn_exp) {
	Trig *t = (Trig *) malloc(sizeof(Trig));
	checkAllocation(t);
	t->trig_fn = trig_fn;
	t->x_coef = x_coef;
	t->x_exp = x_exp;
	t->fn_coef = fn_coef;
	t->fn_exp = fn_exp;
	return t;
}

Inv_Trig *create_inv_trig(int trig_fn, double x_coef, double x_exp, double fn_coef, double fn_exp) {
	Inv_Trig *it = (Inv_Trig *) malloc(sizeof(Inv_Trig));
	checkAllocation(it);
	it->trig_fn = trig_fn;
	it->x_coef = x_coef;
	it->x_exp = x_exp;
	it->fn_coef = fn_coef;
	it->fn_exp = fn_exp;
	return it;
}

Function *create_function(Poly **poly_arr, int poly_count, Exp **exp_arr, int exp_count, Log **log_arr,
	int log_count, Trig **trig_arr, int trig_count, Inv_Trig **inv_trig_arr, int inv_trig_count) {
	Function *f = (Function *) malloc(sizeof(Function));
	checkAllocation(f);
	f->poly_arr = poly_arr;
	f->poly_count = poly_count;
	f->exp_arr = exp_arr;
	f->exp_count = exp_count;
	f->log_arr = log_arr;
	f->log_count = log_count;
	f->trig_arr = trig_arr;
	f->trig_count = trig_count;
	f->inv_trig_arr = inv_trig_arr;
	f->inv_trig_count = inv_trig_count;
	return f;
}

double **create_matrix(int n, int m) {
	int i;
	double **matrix = (double **) malloc(n * sizeof(double *));
	checkAllocation(matrix);

	for (i = 0; i < n; i++) {
		matrix[i] = (double *) malloc(m * sizeof(double));
		checkAllocation(matrix[i]);
	}

	return matrix;
}

/* INPUT FUNCTIONS */
void input_x(double *x_coef, double *x_exp) {
	printf("x's cofactor (x_coef):\n");
	scanf(" %lf", x_coef);
	printf("x's exponent (x_exp):\n");
	scanf(" %lf", x_exp);
}

void input_fn(double *fn_coef, double *fn_exp) {
	printf("Function cofactor (fn_coef):\n");
	scanf(" %lf", fn_coef);
	printf("Function exponent (fn_exp):\n");
	scanf(" %lf", fn_exp);
}

void input_poly(Poly **p) {
	double x_coef, x_exp;

	printf("\nPolynomial: x_coef * x ^ x_exp\n\n");

	input_x(&x_coef, &x_exp);
	*p = create_poly(x_coef, x_exp);
}

void input_exp(Exp **e) {
	double x_coef, x_exp, fn_coef, fn_exp, base;

	printf("\nExponential: fn_coef * (base ^ (x_coef * x ^ x_exp)) ^ fn_exp\n\n");

	input_x(&x_coef, &x_exp);
	input_fn(&fn_coef, &fn_exp);

	printf("Base (base):\n");
	scanf(" %lf", &base);

	*e = create_exp(x_coef, x_exp, fn_coef, fn_exp, base);
}

void input_log(Log **l) {
	double x_coef, x_exp, fn_coef, fn_exp, base;

	printf("\nLogarithmic: fn_coef * (log _ base (x_coef * x ^ x_exp)) ^ fn_exp\n\n");

	input_x(&x_coef, &x_exp);
	input_fn(&fn_coef, &fn_exp);

	printf("Base (base):\n");
	scanf(" %lf", &base);

	*l = create_log(x_coef, x_exp, fn_coef, fn_exp, base);
}

void input_trig(Trig **t) {
	int trig_fn;
	double x_coef, x_exp, fn_coef, fn_exp;

	printf("\nTrigonometric: fn_coef * <trig_fn>(x_coef * x ^ x_exp) ^ fn_exp\n\n");

	do {
		printf("Trigonometric function (trig_fn):\nsin: 0, cos: 1, tan: 2, cot: 3\n");
		scanf(" %d", &trig_fn);
	} while (trig_fn < 0 || trig_fn > 3);

	input_x(&x_coef, &x_exp);
	input_fn(&fn_coef, &fn_exp);

	*t = create_trig(trig_fn, x_coef, x_exp, fn_coef, fn_exp);
}

void input_inv_trig(Inv_Trig **it) {
	int trig_fn;
	double x_coef, x_exp, fn_coef, fn_exp;

	printf("\nInverse trigonometric: fn_coef * arc<trig_fn>(x_coef * x ^ x_exp) ^ fn_exp\n\n");

	do {
		printf("Inverse trigonometric function (trig_fn):\narcsin: 0, arccos: 1, arctan: 2, arccot: 3\n");
		scanf(" %d", &trig_fn);
	} while (trig_fn < 0 || trig_fn > 3);

	input_x(&x_coef, &x_exp);
	input_fn(&fn_coef, &fn_exp);

	*it = create_inv_trig(trig_fn, x_coef, x_exp, fn_coef, fn_exp);
}

void input_function(Function **f) {
	int i;
	int poly_count, exp_count, log_count, trig_count, inv_trig_count;
	Poly **poly_arr;
	Exp **exp_arr;
	Log **log_arr;
	Trig **trig_arr;
	Inv_Trig **inv_trig_arr;
	char *msg;

	/* INPUT POLYNOMIALS */
	printf("Polynomial count:\n");
	scanf(" %d", &poly_count);
	poly_arr = (Poly **) malloc(poly_count * sizeof(Poly *));
	checkAllocation(poly_arr);
	for (i = 0; i < poly_count; i++) {
		input_poly(&poly_arr[i]);
		msg = to_string_poly(poly_arr[i]);
		printf("Added: %s\n", msg);
		free(msg);
	}

	/* INPUT EXPONENTIALS */
	printf("\nExponential count:\n");
	scanf(" %d", &exp_count);
	exp_arr = (Exp **) malloc(exp_count * sizeof(Exp *));
	checkAllocation(exp_arr);
	for (i = 0; i < exp_count; i++) {
		input_exp(&exp_arr[i]);
		msg = to_string_exp(exp_arr[i]);
		printf("Added: %s\n", msg);
		free(msg);
	}

	/* INPUT LOGARITHMICS */
	printf("\nLogarithmic count:\n");
	scanf(" %d", &log_count);
	log_arr = (Log **) malloc(log_count * sizeof(Log *));
	checkAllocation(log_arr);
	for (i = 0; i < log_count; i++) {
		input_log(&log_arr[i]);
		msg = to_string_log(log_arr[i]);
		printf("Added: %s\n", msg);
		free(msg);
	}

	/* INPUT TRIGONOMETRICS */
	printf("\nTrigonometric count:\n");
	scanf(" %d", &trig_count);
	trig_arr = (Trig **) malloc(trig_count * sizeof(Trig *));
	checkAllocation(trig_arr);
	for (i = 0; i < trig_count; i++) {
		input_trig(&trig_arr[i]);
		msg = to_string_trig(trig_arr[i]);
		printf("Added: %s\n", msg);
		free(msg);
	}

	/* INPUT INVERSE TRIGONOMETRICS */
	printf("\nInverse trigonometric count:\n");
	scanf(" %d", &inv_trig_count);
	inv_trig_arr = (Inv_Trig **) malloc(inv_trig_count * sizeof(Inv_Trig *));
	checkAllocation(inv_trig_arr);
	for (i = 0; i < inv_trig_count; i++) {
		input_inv_trig(&inv_trig_arr[i]);
		msg = to_string_inv_trig(inv_trig_arr[i]);
		printf("Added: %s\n", msg);
		free(msg);
	}

	/* CREATE A FUNCTION FROM ALL THE EXPRESSIONS */
	*f = create_function(poly_arr, poly_count, exp_arr, exp_count, log_arr,
		log_count, trig_arr, trig_count, inv_trig_arr, inv_trig_count);

	msg = to_string_function(*f);
	printf("\nFunction:%s\n\n", msg);
	free(msg);
}

void input_matrix(double ***matrix, int n, int m) {
	int i, j;
	*matrix = create_matrix(n, m);
	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			printf("[%d][%d]:\n", i, j);
			scanf(" %lf", &(*matrix)[i][j]);
		}
	}
	printf("Matrix: ");
	print_matrix(*matrix, n, m);
}

/* FUNCTION CALCULATORS */
double calculate_poly(Poly *p, double x) {
	return p->x_coef * pow(x, p->x_exp);
}

double calculate_exp(Exp *e, double x) {
	return e->fn_coef * pow(pow(e->base, e->x_coef * pow(x, e->x_exp)), e->fn_exp);
}

double calculate_log(Log *l, double x) {
	return l->fn_coef * pow(log(l->x_coef * pow(x, l->x_exp)) / log(l->base), l->fn_exp);
}

double calculate_trig(Trig *t, double x) {
	double (*trig_fn)(double);
	switch (t->trig_fn) {
	case 0: trig_fn = sin; break;
	case 1: trig_fn = cos; break;
	case 2: trig_fn = tan; break;
	case 3: trig_fn = cot; break;
	}
	return t->fn_coef * pow(trig_fn(t->x_coef * pow(x, t->x_exp)), t->fn_exp);
}

double calculate_inv_trig(Inv_Trig *it, double x) {
	double (*trig_fn)(double);
	switch (it->trig_fn) {
	case 0: trig_fn = asin; break;
	case 1: trig_fn = acos; break;
	case 2: trig_fn = atan; break;
	case 3: trig_fn = acot; break;
	}
	return it->fn_coef * pow(trig_fn(it->x_coef * pow(x, it->x_exp)), it->fn_exp);
}

double calculate(Function *f, double x) {
	int i;
	double result;
	
	result = 0;
	
	for (i = 0; i < f->poly_count; i++) {
		result += calculate_poly(f->poly_arr[i], x);
	}
	
	for (i = 0; i < f->exp_count; i++) {
		result += calculate_exp(f->exp_arr[i], x);
	}
	
	for (i = 0; i < f->log_count; i++) {
		result += calculate_log(f->log_arr[i], x);
	}
	
	for (i = 0; i < f->trig_count; i++) {
		result += calculate_trig(f->trig_arr[i], x);
	}
	
	for (i = 0; i < f->inv_trig_count; i++) {
		result += calculate_inv_trig(f->inv_trig_arr[i], x);
	}

	return result;
}

/* TO STRING FUNCTIONS */
char *to_string_poly(Poly *p) {
	char *str = (char *) malloc(EXPR_MAX_LENGTH * sizeof(char));
	checkAllocation(str);
	sprintf(str, "%lf * x ^ %lf", p->x_coef, p->x_exp);
	return str;
}

char *to_string_exp(Exp *e) {
	char *str = (char *) malloc(EXPR_MAX_LENGTH * sizeof(char));
	checkAllocation(str);
	sprintf(str, "%lf * (%lf ^ (%lf * x ^ %lf)) ^ %lf",
		e->fn_coef, e->base, e->x_coef, e->x_exp, e->fn_exp);
	return str;
}

char *to_string_log(Log *l) {
	char *str = (char *) malloc(EXPR_MAX_LENGTH * sizeof(char));
	checkAllocation(str);
	sprintf(str, "%lf * (log _ %lf (%lf * x ^ %lf)) ^ %lf",
		l->fn_coef, l->base, l->x_coef, l->x_exp, l->fn_exp);
	return str;
}

char *to_string_trig(Trig *t) {
	char *str = (char *) malloc(EXPR_MAX_LENGTH * sizeof(char));
	checkAllocation(str);
	char trig_fn[4];
	switch (t->trig_fn) {
	case 0: strcpy(trig_fn, "sin"); break;
	case 1: strcpy(trig_fn, "cos"); break;
	case 2: strcpy(trig_fn, "tan"); break;
	case 3: strcpy(trig_fn, "cot"); break;
	}
	sprintf(str, "%lf * %s(%lf * x ^ %lf) ^ %lf",
		t->fn_coef, trig_fn, t->x_coef, t->x_exp, t->fn_exp);
	return str;
}

char *to_string_inv_trig(Inv_Trig *it) {
	char *str = (char *) malloc(EXPR_MAX_LENGTH * sizeof(char));
	checkAllocation(str);
	char trig_fn[7];
	switch (it->trig_fn) {
	case 0: strcpy(trig_fn, "arcsin"); break;
	case 1: strcpy(trig_fn, "arccos"); break;
	case 2: strcpy(trig_fn, "arctan"); break;
	case 3: strcpy(trig_fn, "arccot"); break;
	}
	sprintf(str, "%lf * %s(%lf * x ^ %lf) ^ %lf",
		it->fn_coef, trig_fn, it->x_coef, it->x_exp, it->fn_exp);
	return str;
}

char *to_string_function(Function *f) {
	int i;
	char *str, *strIndex, *add;
	
	str = (char *) malloc((EXPR_MAX_LENGTH * (f->poly_count + f->exp_count +
		f->log_count + f->trig_count + f->inv_trig_count) + 1) * sizeof(char));
	checkAllocation(str);
	strIndex = str;

	for (i = 0; i < f->poly_count; i++) {
		add = to_string_poly(f->poly_arr[i]);
		strIndex += sprintf(strIndex, " %s ", add);
		free(add);
		if (i != f->poly_count - 1 || f->exp_count != 0 || f->log_count != 0 || f->trig_count != 0 || f->inv_trig_count != 0)
			strIndex += sprintf(strIndex, "+");
	}

	for (i = 0; i < f->exp_count; i++) {
		add = to_string_exp(f->exp_arr[i]);
		strIndex += sprintf(strIndex, " %s ", add);
		free(add);
		if (i != f->exp_count - 1 || f->log_count != 0 || f->trig_count != 0 || f->inv_trig_count != 0)
			strIndex += sprintf(strIndex, "+");
	}

	for (i = 0; i < f->log_count; i++) {
		add = to_string_log(f->log_arr[i]);
		strIndex += sprintf(strIndex, " %s ", add);
		free(add);
		if (i != f->log_count - 1 || f->trig_count != 0 || f->inv_trig_count != 0) strIndex += sprintf(strIndex, "+");
	}

	for (i = 0; i < f->trig_count; i++) {
		add = to_string_trig(f->trig_arr[i]);
		strIndex += sprintf(strIndex, " %s ", add);
		free(add);
		if (i != f->trig_count - 1 || f->inv_trig_count != 0) strIndex += sprintf(strIndex, "+");
	}

	for (i = 0; i < f->inv_trig_count; i++) {
		add = to_string_inv_trig(f->inv_trig_arr[i]);
		strIndex += sprintf(strIndex, " %s ", add);
		free(add);
		if (i != f->inv_trig_count - 1) strIndex += sprintf(strIndex, "+");
	}

	return str;
}

/* DESTRUCTORS */
void free_function(Function *f) {
	int i;

	for (i = 0; i < f->poly_count; i++) {
		free(f->poly_arr[i]);
	}
	free(f->poly_arr);

	for (i = 0; i < f->exp_count; i++) {
		free(f->exp_arr[i]);
	}
	free(f->exp_arr);

	for (i = 0; i < f->log_count; i++) {
		free(f->log_arr[i]);
	}
	free(f->log_arr);

	for (i = 0; i < f->trig_count; i++) {
		free(f->trig_arr[i]);
	}
	free(f->trig_arr);

	for (i = 0; i < f->inv_trig_count; i++) {
		free(f->inv_trig_arr[i]);
	}
	free(f->inv_trig_arr);
	free(f);
}

void free_matrix(double **matrix, int n) {
	int i;
	for (i = 0; i < n; i++) {
		free(matrix[i]);
	}
	free(matrix);
}

/* PROCESSES */
double bisection(Function *f, double start, double end, double e, int opt, int max_iter) {
	double mid, fmid, fend;
	int iter;
	
	iter = 0;
	fmid = e + 1;
	
	while (iter < max_iter && ((opt == 1 && fabs(fmid) > e) || (opt == 2 && ((end - start) / pow(2, iter)) > e))) {
		mid = (start + end) / 2;
		fmid = calculate(f, mid);
		fend = calculate(f, end);

		iter++;

		printf("\n%-10s: %+lf\n%-10s: %+lf\n%-10s: %+lf\n%-10s: %+lf\n%-10s: %+lf\n%-10s: %+lf\n%-10s: %d\n",
			"start", start, "end", end, "mid", mid, "f(start)", calculate(f, start), "f(end)", fend, "f(mid)", fmid, "iteration", iter);

		if (fmid * fend <= 0) {
			start = mid;
		} else {
			end = mid;
		}
	}
	
	if (iter == max_iter) {
		printf("Reached max iterations.\n");
	}
	
	return mid;
}

double regula_falsi(Function *f, double start, double end, double e, int opt, int max_iter) {
	double point;
	double fstart, fend, fpoint;
	int iter;
	
	iter = 0;
	fpoint = e + 1;

	while (iter < max_iter && ((opt == 1 && fabs(fpoint) > e) || (opt == 2 && ((end - start) / pow(2, iter)) > e))) {
		fstart = calculate(f, start);
		fend = calculate(f, end);
		point = (end * fstart - start * fend) / (fstart - fend);
		fpoint = calculate(f, point);

		iter++;

		printf("\n%-10s: %+lf\n%-10s: %+lf\n%-10s: %+lf\n%-10s: %+lf\n%-10s: %+lf\n%-10s: %+lf\n%-10s: %d\n",
			"start", start, "end", end, "point", point, "f(start)", fstart, "f(end)", fend, "f(point)", fpoint, "iteration", iter);
		
		if (fstart * fpoint <= 0) {
			end = point;
		} else {
			start = point;
		}
	}
	
	if (iter == max_iter) {
		printf("Reached max iterations.\n");
	}

	return point;
}

double newton_raphson(Function *f, double x0, double e, int max_iter) {
	double xn, xn1, fxn, fpxn;
	int iter;

	xn = xn1 = x0;
	iter = 0;

	do {
		xn = xn1;
		fxn = calculate(f, xn);
		fpxn = centered_difference(f, xn, 0.000001);
		xn1 = xn - fxn / fpxn;

		iter++;

		printf("\n%-10s: %+lf\n%-10s: %+lf\n%-10s: %+lf\n%-10s: %+lf\n%-10s: %d\n",
			"xn", xn, "xn+1", xn1, "f(xn)", fxn, "f'(xn)", fpxn, "iteration", iter);		
	} while (iter < max_iter && fabs(xn1 - xn) > e);
	
	if (iter == max_iter) {
		printf("Reached max iterations.\n");
	}
	
	return xn1;
}

double **inverse_matrix(double **matrix, int n) {
	double **adj = adjoint_matrix(matrix, n);
	double det = determinant_matrix(matrix, n);
	int i, j;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			adj[i][j] /= det;
		}
	}

	return adj;
}

double **gauss_elimination(double **a, double **c, int n) {
	int i, j, k;
	double **x = create_matrix(n, 1);
	double **aug = augmented_matrix(a, c, n);
	double temp;

	for (i = 0; i < n; i++) {
		temp = aug[i][i];
		for (j = i; j < n + 1; j++) {
			aug[i][j] /= temp;
		}

		for (j = i + 1; j < n; j++) {
			temp = aug[j][i];
			for (k = i; k < n + 1; k++) {
				aug[j][k] -= aug[i][k] * temp;
			}
		}
	}

	for (i = n - 1; i >= 0; i--) {
		temp = aug[i][n];
		for (j = n - 1; j > i; j--) {
			temp -= x[j][0] * aug[i][j];
		}
		x[i][0] = temp;
	}

	free_matrix(aug, n);
	return x;
}

double **gauss_seidel(double **a, double **c, double **starting_values, int n, double epsilon, int opt, int max_iter) {
	int i, j;
	double **x = create_matrix(n, 1);
	double **aug = augmented_matrix(a, c, n);
	double **dominant = dominant_matrix(aug, n, n + 1, opt);
	double *prev = (double *) malloc(n * sizeof(double));
	double temp;
	int iter = 0;
	int stop = 0;

	for (i = 0; i < n; i++) {
		x[i][0] = starting_values[i][0];
	}

	while (iter < max_iter && !stop) {
		for (i = 0; i < n; i++) {
			prev[i] = x[i][0];
		}

		for (i = 0; i < n; i++) {
			temp = dominant[i][n];
			for (j = 0; j < n; j++) {
				if (i != j) {
					temp -= dominant[i][j] * x[j][0];
				}
			}
			x[i][0] = temp / dominant[i][i];
		}

		iter++;

		printf("\nvalues: ");
		print_matrix(x, n, 1);
		printf("iteration: %d\n", iter);

		i = 0;
		while (i < n && !stop) {
			if (fabs(x[i][0] - prev[i]) > epsilon) {
				stop = 1;
			} else {
				i++;
			}
		}
		stop = i >= n;
	}
	
	if (iter == max_iter) {
		printf("Reached max iterations.\n");
	}

	free_matrix(aug, n);
	free_matrix(dominant, n);
	free(prev);
	return x;
}

double backward_difference(Function *f, double x, double h) {
	return (calculate(f, x) - calculate(f, x - h)) / h;
}

double forward_difference(Function *f, double x, double h) {
	return (calculate(f, x + h) - calculate(f, x)) / h;
}

double centered_difference(Function *f, double x, double h) {
	return (calculate(f, x + h) - calculate(f, x - h)) / (2 * h);
}

double simpsons_rule_1_3(Function *f, double start, double end, int n) {
	double h = (end - start) / n;
	double result = calculate(f, start) + calculate(f, end);
	int i;

	for (i = 1; i < n; i += 2) {
		result += 4 * calculate(f, start + i * h);
	}
	for (i = 2; i < n - 1; i += 2) {
		result += 2 * calculate(f, start + i * h);
	}

	result *= h / 3;
	
	return result;
}

double simpsons_rule_3_8(Function *f, double start, double end, int n) {
	double h, x1, x2;

	if (n == 1) {
		h = (end - start) / 3;
		x1 = start + h;
		x2 = start + 2 * h;
		
		return (end - start) * (calculate(f, start) + 3 * calculate(f, x1) + 3 * calculate(f, x2) + calculate(f, end)) / 8;
	} else if (n > 1) {
		return simpsons_rule_3_8(f, start, (start + end) / 2, n - 1) + simpsons_rule_3_8(f, (start + end) / 2, end, n - 1);
	} else {
		return 0;
	}
}

double trapezoidal_rule(Function *f, double start, double end, int n) {
	double h = (end - start) / n;
	double result = (calculate(f, start) + calculate(f, end)) / 2;
	int i;

	for (i = 1; i < n; i++) {
		result += calculate(f, start + i * h);
	}

	result *= h;
	
	return result;
}

Function *gregory_newton(int n, double x0, double h, double *vals) {
	int i, j, done = 0, eq;
	double *diffs = (double *) malloc(n * sizeof(double));
	Poly **result, **poly_arr, **acc, **temp, **temp1, **temp2;
	int poly_count;
	double coef;
	Poly *swap;
	
	i = n;
	while (i > 1 && !done) {
		diffs[n - i] = vals[0];

		j = 0;
		eq = 1;
		while (j < i - 1 && eq) {
			if (vals[j] != vals[j + 1]) {
				eq = 0;
			} else {
				j++;
			}
		}

		if (eq) {
			done = 1;
		} else {
			for (j = 0; j < i - 1; j++) {
				vals[j] = vals[j + 1] - vals[j];
			}
			i--;
		}
	}
	diffs[n - 1] = vals[0];
	
	poly_count = n - i + 1;
	result = (Poly **) malloc(1 * sizeof(Poly *));
	result[0] = create_poly(diffs[0], 0);
	acc = (Poly **) malloc(1 * sizeof(Poly *));;
	acc[0] = create_poly(1, 0);

	for (i = 1; i < poly_count; i++) {
		poly_arr = (Poly **) malloc(2 * sizeof(Poly *));
		poly_arr[0] = create_poly(1, 1);
		poly_arr[1] = create_poly(-1 * x0, 0);
		temp = acc;
		acc = multiply_poly_function(poly_arr, 2, acc, i);
		coef = diffs[i] / pow(h, i) / fact(i);
		for (j = 0; j < i + 1; j++) {
			acc[j]->x_coef *= coef;
		}
		temp1 = result;
		result = merge_poly(result, i, acc, i + 1);
		for (j = 0; j < i + 1; j++) {
			acc[j]->x_coef /= coef;
		}
		free(poly_arr[0]);
		free(poly_arr[1]);
		free(poly_arr);
		for (j = 0; j < i; j++) {
			free(temp[j]);
			free(temp1[j]);
		}
		free(temp);
		free(temp1);
		x0 += h;
	}

	for (i = 0; i < poly_count / 2; i++) {
		swap = result[i];
		result[i] = result[poly_count - i - 1];
		result[poly_count - i - 1] = swap;
	}

	free(diffs);
	for (i = 0; i < poly_count; i++) {
		free(acc[i]);
	}
	free(acc);
	return create_function(result, poly_count, NULL, 0, NULL, 0, NULL, 0, NULL, 0);
}

/* HELPER FUNCTIONS */
double cot(double x) {
	return 1.0 / tan(x);
}

double acot(double x) {
	return atan(1.0 / x);
}

void checkAllocation(void *ptr) {
	if (ptr == NULL) {
		fprintf(stderr, "Allocation failed.");
		exit(1);
	}
}

void print_matrix(double **matrix, int n, int m) {
	int i, j;
	printf("[\n");
	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			printf("%12.6lf", matrix[i][j]);
		}
		printf("\n");
	}
	printf("]\n");
}

double determinant_matrix(double **matrix, int n) {
	if (n == 0) return 0;
	if (n == 1) return matrix[0][0];
	if (n == 2) return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];

	int i;
	double sum = 0;
	double **minor;

	for (i = 0; i < n; i++) {
		minor = minor_matrix(matrix, n, 0, i);
		sum += matrix[0][i] * determinant_matrix(minor, n - 1) * pow(-1, i);
		free_matrix(minor, n - 1);
	}

	return sum;
}

double **minor_matrix(double **matrix, int n, int i, int j) {
	int k, l;
	double result;
	double **minor = create_matrix(n - 1, n - 1);

	for (k = 0; k < n; k++) {
		for (l = 0; l < n; l++) {
			if (k != i && l != j) {
				minor[k > i ? k - 1 : k][l > j ? l - 1 : l] = matrix[k][l];
			}
		}
	}

	return minor;
}

double **adjoint_matrix(double **matrix, int n) {
	int i, j;
	double **result = create_matrix(n, n);
	double **minor;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			minor = minor_matrix(matrix, n, i, j);
			result[j][i] = determinant_matrix(minor, n - 1) * pow(-1, i + j);
			free_matrix(minor, n - 1);
		}
	}

	return result;
}

double **augmented_matrix(double **a, double **c, int n) {
	int i, j;
	double **aug = create_matrix(n, n + 1);

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			aug[i][j] = a[i][j];
		}
	}

	for (i = 0; i < n; i++) {
		aug[i][n] = c[i][0];
	}

	return aug;
}

int dominant_matrix_max_values(double **matrix, int n, int row, int col, double current, double *max, double *vals) {
	double a, b;
	double **minor;
	if (col < n) {
		if (n == 2) {
			a = fabs(matrix[0][0] * matrix[1][1]);
			b = fabs(matrix[0][1] * matrix[1][0]);
			if (a >= b) {
				current *= a;
				if (current > *max) {
					vals[row] = matrix[0][0];
					vals[row + 1] = matrix[1][1];
					*max = current;
					return 1;
				}
			} else {
				current *= b;
				if (current > *max) {
					vals[row] = matrix[0][1];
					vals[row + 1] = matrix[1][0];
					*max = current;
					return 1;
				}
			}
			return 0;
		} else if (n > 2) {
			minor = minor_matrix(matrix, n, 0, col);
			a = dominant_matrix_max_values(minor, n - 1, row + 1, 0, fabs(current * matrix[0][col]), max, vals);
			if (a) {
				vals[row] = matrix[0][col];
			}
			free_matrix(minor, n - 1);
			return dominant_matrix_max_values(matrix, n, row, col + 1, current, max, vals) || a;
		}
	}
	return 0;
}

double **dominant_matrix(double **matrix, int n, int m, int opt) {
	int i, j, k, found;
	double **dominant = create_matrix(n, m);
	double *vals = (double *) calloc(n, sizeof(double));
	double *set = (double *) calloc(n, sizeof(double));
	double max = -1;

	if (opt) {
		dominant_matrix_max_values(matrix, n, 0, 0, 1, &max, vals);

		for (i = 0; i < n; i++) {
			found = 0;
			j = 0;
			while (j < n && !found) {
				if (!set[j] && matrix[i][j] == vals[i]) {
					found = 1;
				} else {
					j++;
				}
			}
			if (found) {
				for (k = 0; k < m; k++) {
					dominant[j][k] = matrix[i][k];
				}
				set[j] = 1;
			}
		}
	}
	
	i = 0;
	while (i < n && set[i]) {
		i++;
	}

	if (i < n) {
		if (opt) printf("Couldn't find diagonally dominant matrix.\n");
		for (i = 0; i < n; i++) {
			for (j = 0; j < m; j++) {
				dominant[i][j] = matrix[i][j];
			}
		}
	} else {
		printf("Diagonally dominant matrix: ");
		print_matrix(dominant, n, n);
	}
	
	free(vals);
	free(set);
	return dominant;
}

Poly *multiply_poly(Poly *p1, Poly *p2) {
	return create_poly(p1->x_coef * p2->x_coef, p1->x_exp + p2->x_exp);
}

Poly **multiply_poly_function(Poly **poly_arr1, int n1, Poly **poly_arr2, int n2) {
	int n3 = n1 + n2 - 1;
	Poly **result = (Poly **) malloc(n3 * sizeof(Poly *));
	Poly *temp;
	int i, j, k;

	checkAllocation(result);
	for (i = 0; i < n3; i++) {
		result[i] = NULL;
	}

	for (i = 0; i < n1; i++) {
		if (poly_arr1[i] != NULL) {
			for (j = 0; j < n2; j++) {
				if (poly_arr2[j] != NULL) {
					temp = multiply_poly(poly_arr1[i], poly_arr2[j]);
					k = 0;
					while (k < n3 && result[k] != NULL && result[k]->x_exp != temp->x_exp) {
						k++;
					}
					if (k < n3) {
						if (result[k] == NULL) {
							result[k] = temp;
						} else if (result[k]->x_exp == temp->x_exp) {
							result[k]->x_coef += temp->x_coef;
							free(temp);
						}
					}
				}
			}
		}
	}

	return result;
}

Poly **merge_poly(Poly **poly_arr1, int n1, Poly **poly_arr2, int n2) {
	Poly **simplified;
	int i, j, k, eq = 0, found;

	for (i = 0; i < n1; i++) {
		for (j = 0; j < n2; j++) {
			if (poly_arr1[i] != NULL && poly_arr2[j] != NULL && poly_arr1[i]->x_exp == poly_arr2[j]->x_exp) {
				eq++;
			}
		}
	}

	simplified = (Poly **) malloc((n1 + n2 - eq) * sizeof(Poly *));
	k = 0;

	for (i = 0; i < n1; i++) {
		if (poly_arr1[i] != NULL) {
			simplified[k++] = create_poly(poly_arr1[i]->x_coef, poly_arr1[i]->x_exp);
		}
	}

	for (i = 0; i < n2; i++) {
		if (poly_arr2[i] != NULL) {
			found = 0;
			j = 0;
			while (j < n1 && !found) {
				if (poly_arr1[j] != NULL && poly_arr2[i]->x_exp == poly_arr1[j]->x_exp) {
					found = 1;
				} else {
					j++;
				}
			}
			if (found) {
				simplified[j]->x_coef += poly_arr2[i]->x_coef;
			} else {
				simplified[k++] = create_poly(poly_arr2[i]->x_coef, poly_arr2[i]->x_exp);
			}
		}
	}

	return simplified;
}

double fact(int n) {
	if (n < 0) return 0;
	if (n < 2) return 1;
	int result = n;
	while (--n > 1) result *= n;
	return result;
}
