#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define PI 3.14159265358979323846

void fft_recursive(complex double *arr, int n, int invert) {
    if (n == 1) return;

    complex double *even = malloc(n / 2 * sizeof(complex double));
    complex double *odd = malloc(n / 2 * sizeof(complex double));

    for (int i = 0; i < n / 2; i++) {
        even[i] = arr[i * 2];
        odd[i] = arr[i * 2 + 1];
    }

    fft_recursive(even, n / 2, invert);
    fft_recursive(odd, n / 2, invert);

    double angle = 2 * PI / n * (invert ? -1 : 1);
    complex double wn = cexp(angle * I);
    complex double w = 1;

    for (int i = 0; i < n / 2; i++) {
        arr[i] = even[i] + w * odd[i];
        arr[i + n / 2] = even[i] - w * odd[i];
        if (invert) {
            arr[i] /= 2;
            arr[i + n / 2] /= 2;
        }
        w *= wn;
    }

    free(even);
    free(odd);
}

void fft_convolution(int *x, int *y, int n, int *result) {
    complex double *x_fft = malloc(n * sizeof(complex double));
    complex double *y_fft = malloc(n * sizeof(complex double));

    for (int i = 0; i < n; i++) {
        x_fft[i] = i < n ? x[i] : 0;
        y_fft[i] = i < n ? y[i] : 0;
    }

    fft_recursive(x_fft, n, 0);
    fft_recursive(y_fft, n, 0);

    complex double *result_fft = malloc(n * sizeof(complex double));
    for (int i = 0; i < n; i++) {
        result_fft[i] = x_fft[i] * y_fft[i];
    }

    fft_recursive(result_fft, n, 1);

    for (int i = 0; i < n; i++) {
        result[i] = (int)round(creal(result_fft[i]));
    }

    // Handle carry-over
    for (int i = 0; i < n - 1; i++) {
        if (result[i] >= 10) {
            result[i + 1] += result[i] / 10;
            result[i] %= 10;
        }
    }

    free(x_fft);
    free(y_fft);
    free(result_fft);
}

void to_coeff_array(int num, int *coeff_array, int *size) {
    *size = 0;
    while (num > 0) {
        coeff_array[(*size)++] = num % 10;
        num /= 10;
    }
}

int from_coeff_array(int *coeff_array, int size) {
    int result = 0;
    for (int i = size - 1; i >= 0; i--) {
        result = result * 10 + coeff_array[i];
    }
    return result;
}

int fft_multiply(int a, int b) {
    int coeff_a[1024] = {0};
    int coeff_b[1024] = {0};
    int size_a, size_b;

    to_coeff_array(a, coeff_a, &size_a);
    to_coeff_array(b, coeff_b, &size_b);

    int n = 1;
    while (n < size_a + size_b - 1) {
        n <<= 1;
    }

    int *product_coeff = malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) product_coeff[i] = 0;

    fft_convolution(coeff_a, coeff_b, n, product_coeff);

    int result = from_coeff_array(product_coeff, n);

    if ((a < 0) ^ (b < 0)) {
        result = -result;
    }

    free(product_coeff);

    return result;
}

int main() {
    int x = 12345;
    int y = 6789;
    printf("%d * %d = %d\n", x, y, fft_multiply(x, y));
    return 0;
}
