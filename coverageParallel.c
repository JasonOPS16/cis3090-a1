#define _POSIX_C_SOURCE 200809L // to fix the implicit declaration warning

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef struct {
    double width, height;
    long k; // num nodes
    double *x; // x pos
    double *y; // y pos
    double *r2; // r^2
} Field;

// parse for unsigned long long, rejecting anything above max, 0 for success -1 for error
static int parseUnsigned(const char *s, unsigned long long max, unsigned long long *out) {
    if (*s == '\0') {
        return -1; // empty string
    }
    unsigned long long v = 0;
    for (; *s != '\0'; s++) {
        if (*s < '0' || *s > '9') {
            return -1; // not a digit
        }
        unsigned long long d = (unsigned long long)(*s - '0');
        if (v > (max - d) / 10) {
            return -1; // exceeds max
        }
        v = v * 10 + d;
    }
    *out = v;
    return 0; 
}

// free all mem for field and set to NULL
static void freeField(Field *f) {
    free(f->x);
    free(f->y);
    free(f->r2);
    f->x = f->y = f->r2 = NULL;
}

// load field from file 0 or pass -1 on error
static int loadField(const char *path, Field *f) {
    f->x = f->y = f->r2 = NULL;

    FILE *fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "error cant open input file '%s'\n", path);
        return -1;
    }

    // read width, height
    if (fscanf(fp, "%lf %lf", &f->width, &f->height) != 2 || f->width <= 0 || f->height <= 0) {
        fprintf(stderr, "error bad field size in '%s'\n", path);
        fclose(fp);
        return -1;
    }

    // read k
    if (fscanf(fp, "%ld", &f->k) != 1 || f->k < 0) {
        fprintf(stderr, "error bade node count in '%s'\n", path);
        fclose(fp);
        return -1;
    }

    // allocate arrays for k nodes
    if (f->k > 0) {
        f->x = malloc((size_t)f->k * sizeof(double));
        f->y = malloc((size_t)f->k * sizeof(double));
        f->r2 = malloc((size_t)f->k * sizeof(double));
        if (!f->x || !f->y || !f->r2) {
            fprintf(stderr, "error outof mem\n");
            freeField(f);
            fclose(fp);
            return -1;
        }
    }

    // read k lines of x y r
    for (long i = 0; i < f->k; i++) {
        double r;
        if (fscanf(fp, "%lf %lf %lf", &f->x[i], &f->y[i], &r) != 3 || r < 0) {
            fprintf(stderr, "error bad node line %ld in '%s'\n", i, path);
            freeField(f);
            fclose(fp);
            return -1;
        }
        f->r2[i] = r * r;
    }

    fclose(fp);
    return 0;
}

// uniform double in [0,1)
static inline double uniform01(unsigned int *seed) {
    return (double)rand_r(seed) / ((double)RAND_MAX + 1.0);
}

// covered if it is in a rad of at least one node
static inline int isCovered(const Field *f, double px, double py) {
    for (long i = 0; i < f->k; i++) {
        double dx = px - f->x[i];
        double dy = py - f->y[i];
        if (dx * dx + dy * dy < f->r2[i]) return 1; /* stop at first hit */
    }
    return 0;
}

// return elapsed seconds
static double elapsedSeconds(const struct timespec *a, const struct timespec *b) {
    return (double)(b->tv_sec - a->tv_sec) + (double)(b->tv_nsec - a->tv_nsec) / 1e9;
}

// using the foster's methodology:
// partition: sample points are independent, so each thread throws its part
// serial is only to add up in main

// communicate: each thread gets its count, seed and the read only field
// in its own task struct

// agglomerate: each thread gets N/p points, one per thread is too small

// map: thread i gets n/p points, +1 for the first n%p threads thefore total n,
// each thread gets its own seed baseSeed + i


// task struct for each thread
typedef struct {
    const Field *field; // read only and shared
    long long samples; // num of points to throw
    unsigned int seed; // baseSeed + i
    long long hits; // out put
} Task;

// throws samples and counts the hits
static void *worker(void *arg) {
    Task *t = (Task *)arg;

    // local counter
    long long hits = 0;
    for (long long s = 0; s < t->samples; s++) {
        double px = uniform01(&t->seed) * t->field->width;
        double py = uniform01(&t->seed) * t->field->height;
        hits += isCovered(t->field, px, py);
    }

    // write to threads task stuct main will read after p join
    t->hits = hits;

    return NULL;
}

int main(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <inputFile> <N> <p> <baseSeed>\n", argv[0]);
        return 1;
    }

    unsigned long long nParsed, pParsed, seedParsed;
    if (parseUnsigned(argv[2], 1000000000000000000ULL, &nParsed) != 0 || nParsed == 0) {
        fprintf(stderr, "error n must be a pos int <= max 10^18\n");
        return 1;
    }

    // p must be at least 1, cap at 64
    if (parseUnsigned(argv[3], 64ULL, &pParsed) != 0 || pParsed == 0) {
        fprintf(stderr, "error p must be a pos int <= max 64\n");
        return 1;
    }

    // check if fits in unsigned int range
    if (parseUnsigned(argv[4], (unsigned int)-1, &seedParsed) != 0) {
        fprintf(stderr, "error baseSeed has to be non neg int that fits a unsigned int\n");
        return 1;
    }
    long long n = (long long)nParsed;
    int p = (int)pParsed;
    unsigned int baseSeed = (unsigned int)seedParsed;

    // load field from file
    Field field;
    if (loadField(argv[1], &field) != 0) return 1;

    // alloc threads and tasks
    pthread_t *threads = malloc((size_t)p * sizeof(pthread_t));
    Task *tasks = malloc((size_t)p * sizeof(Task));
    if (!threads || !tasks) {
        fprintf(stderr, "error out of mem\n");
        free(threads);
        free(tasks);
        freeField(&field);
        return 1;
    }

    // fill tasks struct for each thread the mapping step
    long long base = n / p;
    long long extra = n % p; // first 'extra' threads get one more point
    for (int i = 0; i < p; i++) {
        tasks[i].field = &field;
        tasks[i].samples = base + (i < extra ? 1 : 0);
        tasks[i].seed = baseSeed + (unsigned int)i;
        tasks[i].hits = 0;
    }

    // timed area
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    // create the threads
    int created = 0;
    for (int i = 0; i < p; i++) {
        if (pthread_create(&threads[i], NULL, worker, &tasks[i]) != 0) {
            fprintf(stderr, "error cant create thread %d\n", i);
            break;
        }
        created++;
    }

    // wait for all threads that started
    for (int i = 0; i < created; i++) {
        pthread_join(threads[i], NULL);
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    // end of timed area

    // thread creation failed failled
    if (created != p) {
        free(threads);
        free(tasks);
        freeField(&field);
        return 1;
    }

    // serial combine step
    long long hits = 0;
    for (int i = 0; i < p; i++) {
        hits += tasks[i].hits;
    }

    double fraction = (double)hits / (double)n;
    double area = fraction * field.width * field.height;

    printf("SAMPLES: %lld\n", n);
    printf("HITS: %lld\n", hits);
    printf("COVERAGE: %.1f %.6f\n", area, fraction * 100.0);
    fprintf(stderr, "TIME: %.6f\n", elapsedSeconds(&start, &end));

    free(threads);
    free(tasks);
    freeField(&field);
    return 0;
}