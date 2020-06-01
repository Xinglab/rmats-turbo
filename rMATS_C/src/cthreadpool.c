#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include "../include/type.h"
#include "../include/global.h"
#include "../include/cthreadpool.h"


void *pool_thread(void *arg) {
    threadpool_t *pool = (threadpool_t*)arg;
    task_t task;
    int head;

    while(1) {
        pthread_mutex_lock(&(pool->pool_lock));

        while(pool->task_count == 0 && pool->flag != POOL_RECLAIM) {
            pthread_cond_wait(&(pool->notify), &(pool->pool_lock));
        }

        if (pool->task_count == 0 && pool->flag == POOL_RECLAIM) {
            break;
        }

        head = pool->task_head;
        task.func = pool->tasks[head].func;
        task.arg = pool->tasks[head].arg;
        pool->task_head++;
        pool->task_count--;

        pthread_mutex_unlock(&(pool->pool_lock));

        *(pool->tasks[head].ret) = (task.func)(task.arg);
    }

    pthread_mutex_unlock(&(pool->pool_lock));

    return NULL;
}


threadpool_t* threadpool_create(const int nthread, const int queue_size) {
    int i = 0;
    threadpool_t *pool = NULL;

    if ((pool = (threadpool_t *)malloc(sizeof(threadpool_t))) == NULL) {
        return NULL;
    }

    pool->nthread = nthread;
    pool->queue_size = queue_size;
    pool->task_head = 0;
    pool->task_tail = -1;
    pool->task_count = 0;

    pthread_mutex_init(&pool->pool_lock, NULL);
    pthread_cond_init(&pool->notify, NULL);

    pool->working = (int*)malloc(sizeof(int)*nthread);
    pool->threads = (pthread_t*)malloc(sizeof(pthread_t)*nthread);
    pool->tasks = (task_t*)malloc(sizeof(task_t)*queue_size);
    for (i = 0; i < nthread; ++i) {
        pool->working[i] = IDLE;
    }

    for(i = 0; i < nthread; i++) {
        if(pthread_create(&(pool->threads[i]), NULL, pool_thread, (void*)pool) != 0) {
            threadpool_reclaim(pool);
            return NULL;
        }
        pool->working[i] = INCOMPLETED;
    }

    return pool;
}


int threadpool_add(threadpool_t *pool, void* (*func)(void *), void *arg, void **ret) {
    int tail;

    if (pool == NULL || func == NULL) {
        return POOL_ADD_ERR;
    }

    if(pthread_mutex_lock(&(pool->pool_lock)) != 0) {
        return POOL_LOCK_FAILURE;
    }

    tail = pool->task_tail + 1;
    tail = tail == pool->queue_size? 0:tail;
    if (pool->task_count == pool->queue_size) {
        // TODO erroneous.
        if ((realloc(pool->tasks, 2 * pool->queue_size)) != 0) {
            return POOL_ADD_ERR;
        }
        tail = pool->queue_size;
        pool->queue_size *= 2;
    }

    pool->tasks[tail].func = func;
    pool->tasks[tail].arg = arg;
    pool->tasks[tail].ret = ret;
    pool->task_tail = tail;
    pool->task_count++;

    if(pthread_cond_signal(&(pool->notify)) != 0 ||
       pthread_mutex_unlock(&pool->pool_lock) != 0) {
        return POOL_LOCK_FAILURE;
    }

    return 0;
}


int threadpool_reclaim(threadpool_t *pool) {
    int i = 0;

    if(pthread_mutex_lock(&(pool->pool_lock)) != 0) {
        return POOL_LOCK_FAILURE;
    }

    pool->flag = POOL_RECLAIM;
    pthread_mutex_unlock(&(pool->pool_lock));

    for (i = 0; i < pool->nthread; ++i) {
        pool->working[i] != IDLE && pthread_join(pool->threads[i], NULL);
    }

    free(pool->threads);
    free(pool->tasks);
    free(pool);

    return 0;
}
