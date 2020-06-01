#ifndef TRANS_CTHREADPOOL
#define TRANS_CTHREADPOOL


#include <pthread.h>
#include "../include/type.h"


typedef struct {
    void *(*func)(void *);
    void *arg;
    void **ret;
} task_t;

typedef struct {
    int nthread;
    int queue_size;
    int task_head;
    int task_tail;
    int task_count;
    int flag;
    int *working;
    pthread_mutex_t pool_lock;
    pthread_cond_t notify;
    pthread_t *threads;
    task_t *tasks;
} threadpool_t;

void update_thread_state(void *idx);

void c_thread_pool(int nthread, pthread_t *ptp_array, void **value_array,
                   void* (*func) (void*), void* arg, int flag);

threadpool_t* threadpool_create(const int nthread, const int queue_size);

int threadpool_add(threadpool_t *pool, void* (*func)(void *), void *arg, void **ret);

void *pool_thread(void *arg);

int threadpool_reclaim(threadpool_t *pool);

#endif
