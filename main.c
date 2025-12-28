#include <stdio.h>
#include <math.h>
#include <locale.h>
#define EPS 1e-8           // точность для численных процедур
#define ZERO_TOL 1e-12     // порог для проверки x=0 (там разрыв)

double f_value(double x, int* ok);
double f_deriv(double x, int* ok);
void print_table(double a, double b, double h);
void min_max(double a, double b, double h);
int bisection(double a, double b, double y, double eps, double* root);
void find_x_for_y(double a, double b, double h, double y, double eps);

int crosses_undefined(double a, double b) {
    /* единственная проблемная точка у функции — x=0 */
    return (a < 0.0 && b > 0.0) || fabs(a) < ZERO_TOL || fabs(b) < ZERO_TOL;
}

/* --- Значение функции f(x). ok=1 если определена, иначе 0 --- */
double f_value(double x, int* ok) {
    *ok = 1;
    if (x < -1.0) {
        return exp(-x * x);
    }
    else if (x >= -1.0 && x < 1.0) {
        if (fabs(x) < ZERO_TOL) { *ok = 0; return 0.0; }
        return log(fabs(x)) / x;
    }
    else { /* x >= 1 */
        double s = 0.0;
        for (int n = 0; n <= 5; ++n) {
            int sign = (n % 2 == 0) ? 1 : -1;           /* (-1)^n */
            s += sign * pow(x, n) / sqrt(n + 2.0);
        }
        return s;
    }
}

/* --- Аналитическая производная f'(x). ok=1 если определена --- */
double f_deriv(double x, int* ok) {
    *ok = 1;
    if (x < -1.0) {
        /* (e^{-x^2})' = -2x e^{-x^2} */
        return -2.0 * x * exp(-x * x);
    }
    else if (x >= -1.0 && x < 1.0) {
        if (fabs(x) < ZERO_TOL) { *ok = 0; return 0.0; }
        /* (ln|x|/x)' = (1 - ln|x|)/x^2 */
        return (1.0 - log(fabs(x))) / (x * x);
    }
    else { /* x >= 1 */
        double s = 0.0;
        for (int n = 1; n <= 5; ++n) {
            int sign = (n % 2 == 0) ? 1 : -1;           /* (-1)^n */
            s += sign * n * pow(x, n - 1) / sqrt(n + 2.0);
        }
        return s;
    }
}

/* --- Таблица x -> f(x) --- */
void print_table(double a, double b, double h) {
    if (h <= 0) { printf("Шаг должен быть > 0.\n"); return; }
    if (a > b) { double t = a; a = b; b = t; }

    printf("\n     x           f(x)\n");
    printf("---------------------------\n");
    for (double x = a; x <= b + 1e-12; x += h) {
        int ok = 0;
        double y = f_value(x, &ok);
        if (ok) printf("%10.6f   %12.6f\n", x, y);
        else    printf("%10.6f   %12s\n", x, "undef");
    }
    printf("---------------------------\n\n");
}

/* --- Глобальные min и max на отрезке (простой перебор по сетке) --- */
void min_max(double a, double b, double h) {
    if (h <= 0) { printf("Шаг должен быть > 0.\n"); return; }
    if (a > b) { double t = a; a = b; b = t; }

    double minv = 1e308, maxv = -1e308;
    double xmin = a, xmax = a;
    int found = 0;

    for (double x = a; x <= b + 1e-12; x += h) {
        int ok = 0;
        double y = f_value(x, &ok);
        if (!ok) continue;                  /* пропускаем x=0 */
        if (!found) { minv = maxv = y; xmin = xmax = x; found = 1; }
        if (y < minv) { minv = y; xmin = x; }
        if (y > maxv) { maxv = y; xmax = x; }
    }

    if (!found) {
        printf("На отрезке нет точек из области определения.\n");
        return;
    }
    printf("\nМинимум: f(%.6f) = %.6f\n", xmin, minv);
    printf("Максимум: f(%.6f) = %.6f\n\n", xmax, maxv);
}

/* --- Бисекция на отрезке [a,b] для уравнения f(x)=y --- */
int bisection(double a, double b, double y, double eps, double* root) {
    int ok1 = 0, ok2 = 0;
    double fa = f_value(a, &ok1) - y;
    double fb = f_value(b, &ok2) - y;

    if (!(ok1 && ok2)) return 0;            /* края не в области определения */
    if (fa * fb > 0)    return 0;           /* нет смены знака */

    for (int it = 0; it < 1000 && fabs(b - a) > eps; ++it) {
        double c = 0.5 * (a + b);
        int okc = 0;
        double fc = f_value(c, &okc) - y;
        if (!okc) {                         /* только если c попало в x=0 */
            c += 1e-6;
            okc = 0; fc = f_value(c, &okc) - y;
            if (!okc) return 0;
        }
        if (fa * fc <= 0) { b = c; fb = fc; }
        else { a = c; fa = fc; }
    }
    *root = 0.5 * (a + b);
    return 1;
}

/* --- Поиск всех решений на [a,b] по сетке: сканируем и сужение бискцией --- */
void find_x_for_y(double a, double b, double h, double y, double eps) {
    if (h <= 0) { printf("Шаг должен быть > 0.\n"); return; }
    if (a > b) { double t = a; a = b; b = t; }

    int found = 0;
    for (double left = a; left < b; left += h) {
        double right = left + h;
        if (right > b) right = b;
        if (crosses_undefined(left, right)) continue;   /* избегаем x=0 */
        int ok1 = 0, ok2 = 0;
        double g1 = f_value(left, &ok1) - y;
        double g2 = f_value(right, &ok2) - y;
        if (!(ok1 && ok2)) continue;

        if (g1 == 0.0) { printf("x = %.8f (точное совпадение)\n", left); found = 1; }
        if (g1 * g2 > 0) continue;

        double r = 0.0;
        if (bisection(left, right, y, eps, &r)) {
            printf("x ~= %.8f  (на отрезке [%.6f, %.6f])\n", r, left, right);
            found = 1;
        }
    }
    if (!found) {
        printf("На данном интервале корней не найдено (для заданного Y).\n");
    }
}

/* ----------------------- MAIN ----------------------- */
int main(void) {
    setlocale(LC_ALL, "");
    while (1) {
        int choice;
        printf("Меню:\n");
        printf("1) Значение f(x) в точке\n");
        printf("2) Таблица x -> f(x) на интервале\n");
        printf("3) Min/Max на отрезке (по сетке)\n");
        printf("4) Поиск x: f(x) ~= Y (бисекция)\n");
        printf("5) Производная f'(x) в точке\n");
        printf("0) Выход\n");
        printf("Ваш выбор: ");
        if (scanf("%d", &choice) != 1) return 0;

        if (choice == 0) break;

        if (choice == 1) {
            double x; printf("Введите x: ");
            scanf("%lf", &x);
            int ok = 0; double y = f_value(x, &ok);
            if (ok) printf("f(%.6f) = %.10f\n\n", x, y);
            else    printf("В точке x=%.6f функция не определена (x=0 в средней ветви).\n\n", x);

        }
        else if (choice == 2) {
            double a, b, h;
            printf("Введите левую границу a, правую b и шаг h: ");
            scanf("%lf %lf %lf", &a, &b, &h);
            print_table(a, b, h);

        }
        else if (choice == 3) {
            double a, b, h;
            printf("Введите левую границу a, правую b и шаг h: ");
            scanf("%lf %lf %lf", &a, &b, &h);
            min_max(a, b, h);

        }
        else if (choice == 4) {
            double a, b, h, y, eps;
            printf("Введите Y, затем a b (интервал поиска), шаг сетки h и точность eps: ");
            scanf("%lf %lf %lf %lf %lf", &y, &a, &b, &h, &eps);
            printf("Ищем решения f(x)=%.10f на [%.6f, %.6f]\n", y, a, b);
            find_x_for_y(a, b, h, y, eps);

        }
        else if (choice == 5) {
            double x; printf("Введите x: ");
            scanf("%lf", &x);
            int ok = 0; double df = f_deriv(x, &ok);
            if (ok) printf("f'(%.6f) = %.10f\n\n", x, df);
            else    printf("Производная не определена в x=%.6f (точка разрыва средней ветви).\n\n", x);

        }
        else {
            printf("Неизвестный пункт меню.\n\n");
        }
    }
    return 0;
}

