/*
This file is part of ``kdtree'', a library for working with kd-trees.
Copyright (C) 2007-2011 John Tsiombikas <nuclear@member.fsf.org>

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. The name of the author may not be used to endorse or promote products
   derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
OF SUCH DAMAGE.
*/
/* single nearest neighbor search written by Tamas Nepusz <tamas@cs.rhul.ac.uk> */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "kdtree.h"

#if defined(WIN32) || defined(__WIN32__)
#include <malloc.h>
#endif

#ifdef USE_LIST_NODE_ALLOCATOR

#ifndef NO_PTHREADS
#include <pthread.h>
#else

#ifndef I_WANT_THREAD_BUGS
#error "You are compiling with the fast list node allocator, with pthreads disabled! This WILL break if used from multiple threads."
#endif	/* I want thread bugs */

#endif	/* pthread support */
#endif	/* use list node allocator */



#define SQ(x)			((x) * (x))
#define HP_LCHILD(x)    (2*(x) + 1)
#define HP_RCHILD(x)    (2*(x) + 2)
#define HP_PARENT(x)    (((x)-1)/2)


static void clear_rec(struct kdnode *node, void (*destr)(void*));
static int insert_rec(struct kdnode **node, const double *pos, void *data, int dir, int dim);
static void clear_results(kdres *set);

static struct kdhyperrect* hyperrect_create(int dim, const double *min, const double *max);
static void hyperrect_free(struct kdhyperrect *rect);
static struct kdhyperrect* hyperrect_duplicate(const struct kdhyperrect *rect);
static void hyperrect_extend(struct kdhyperrect *rect, const double *pos);
static double hyperrect_dist_sq(struct kdhyperrect *rect, const double *pos);
static float dist_sq( float *a1, float *a2, int dims ) ;

static kdres *kd_hp_create(int size);
static void kd_hp_swap(res_node *node1, res_node *node2);
static void kd_heapify(kdres *heap, int i);
static int kd_hp_insert(kdres *heap, res_node node);
static int kd_hp_modify(kdres *heap, res_node node);
static void kd_hp_sort(kdres *heap);
#ifdef USE_LIST_NODE_ALLOCATOR
static struct res_node *alloc_resnode(void);
static void free_resnode(struct res_node*);
#else
#define alloc_resnode()		malloc(sizeof(struct res_node))
#define free_resnode(n)		free(n)
#endif

static float dist_sq( float *a1, float *a2, int dims )
{
  float dist_sq = 0, diff;
  while( --dims >= 0 ) {
    diff = (a1[dims] - a2[dims]);
    dist_sq += diff*diff;
  }
  return dist_sq;
}

struct kdtree *kd_create(int k)
{
    struct kdtree *tree = NULL;

    if(!(tree = (struct kdtree *)malloc(sizeof *tree))) {
        return 0;
    }

    tree->dim = k;
    tree->root = NULL;
    tree->destr = NULL;
    tree->rect = NULL;

    return tree;
}

void kd_free(struct kdtree *tree)
{
	if(tree) {
        kd_clear(tree);
        free(tree);
        tree = NULL;
    }
}

static void clear_rec(struct kdnode *node, void (*destr)(void*))
{
    if(!node)
        return;

    clear_rec(node->left, destr);
    clear_rec(node->right, destr);

    if(destr) {
        destr(node->data);
    }
    free(node->pos);
    node->pos = NULL;
    free(node);
    node = NULL;

}

void kd_clear(struct kdtree *tree)
{
	clear_rec(tree->root, tree->destr);
    tree->root = NULL;
    tree->dim = 0;
    tree->destr = NULL;
	if (tree->rect) {
		hyperrect_free(tree->rect);
        tree->rect = NULL;
	}
}

void kd_data_destructor(struct kdtree *tree, void (*destr)(void*))
{
	tree->destr = destr;
}


static int insert_rec(struct kdnode **nptr, const double *pos, void *data, int dir, int dim)
{
	int new_dir;
	struct kdnode *node;

	if(!*nptr) {
        if(!(node = (struct kdnode *)calloc(1,sizeof *node))) {
			return -1;
		}
        if(!(node->pos = (double*)calloc(dim , sizeof *node->pos))) {
			free(node);
			return -1;
		}
		memcpy(node->pos, pos, dim * sizeof *node->pos);
		node->data = data;
		node->dir = dir;
        node->left = NULL;
        node->right = NULL;
		*nptr = node;
		return 0;
	}

	node = *nptr;
	new_dir = (node->dir + 1) % dim;
	if(pos[node->dir] < node->pos[node->dir]) {
		return insert_rec(&(*nptr)->left, pos, data, new_dir, dim);
	}
	return insert_rec(&(*nptr)->right, pos, data, new_dir, dim);
}

int kd_insert(struct kdtree *tree, const double *pos, void *data)
{
	if (insert_rec(&tree->root, pos, data, 0, tree->dim)) {
		return -1;
	}

	if (tree->rect == 0) {
		tree->rect = hyperrect_create(tree->dim, pos, pos);
	} else {
		hyperrect_extend(tree->rect, pos);
	}

	return 0;
}

int kd_insertf(struct kdtree *tree, const float *pos, void *data)
{
	static double sbuf[16];
	double *bptr, *buf = 0;
	int res, dim = tree->dim;

	if(dim > 16) {
#ifndef NO_ALLOCA
		if(dim <= 256)
			bptr = buf = (double*)alloca(dim * sizeof *bptr);
		else
#endif
			if(!(bptr = buf = (double*)malloc(dim * sizeof *bptr))) {
				return -1;
			}
	} else {
		bptr = buf = sbuf;
	}

	while(dim-- > 0) {
		*bptr++ = *pos++;
	}

	res = kd_insert(tree, buf, data);
#ifndef NO_ALLOCA
	if(tree->dim > 256)
#else
	if(tree->dim > 16)
#endif
		free(buf);
	return res;
}

int kd_insert3(struct kdtree *tree, double x, double y, double z, void *data)
{
	double buf[3];
	buf[0] = x;
	buf[1] = y;
	buf[2] = z;
	return kd_insert(tree, buf, data);
}

int kd_insert3f(struct kdtree *tree, float x, float y, float z, void *data)
{
	double buf[3];
	buf[0] = x;
	buf[1] = y;
	buf[2] = z;
	return kd_insert(tree, buf, data);
}

static void clear_results(kdres *rset)
{
    free(rset->data);
    rset->data = NULL;
    rset->size = 0;

}
void kd_res_free(kdres *rset)
{
    clear_results(rset);
	free(rset);
    rset = NULL;
}

int kd_res_size(kdres *rset)
{
    if(rset)
    {
        return (rset->size);
    }
    else {
        return -1;
    }

}

double kd_res_dis(kdres *rset, int i)
{
    if(rset && i<rset->size)
    {
        return sqrt(rset->data[i].dist_sq);
    }
    else {
        return -1;
    }

}

int kd_res_item(kdres *rset, int i)
{
    if(rset && i < rset->size)
    {
        return (int)rset->data[i].item->data;
    }

}

/* ---- hyperrectangle helpers ---- */
static struct kdhyperrect* hyperrect_create(int dim, const double *min, const double *max)
{
	size_t size = dim * sizeof(double);
	struct kdhyperrect* rect = 0;

    if (!(rect = (struct kdhyperrect*)malloc(sizeof(struct kdhyperrect)))) {
		return 0;
	}

	rect->dim = dim;
    if (!(rect->min = (double*)malloc(size))) {
		free(rect);
		return 0;
	}
    if (!(rect->max = (double*)malloc(size))) {
		free(rect->min);
		free(rect);
		return 0;
	}
	memcpy(rect->min, min, size);
	memcpy(rect->max, max, size);

	return rect;
}

static void hyperrect_free(struct kdhyperrect *rect)
{
	free(rect->min);
	free(rect->max);
	free(rect);
}

static struct kdhyperrect* hyperrect_duplicate(const struct kdhyperrect *rect)
{
	return hyperrect_create(rect->dim, rect->min, rect->max);
}

static void hyperrect_extend(struct kdhyperrect *rect, const double *pos)
{
	int i;

	for (i=0; i < rect->dim; i++) {
		if (pos[i] < rect->min[i]) {
			rect->min[i] = pos[i];
		}
		if (pos[i] > rect->max[i]) {
			rect->max[i] = pos[i];
		}
	}
}

static double hyperrect_dist_sq(struct kdhyperrect *rect, const double *pos)
{
	int i;
	double result = 0;

	for (i=0; i < rect->dim; i++) {
		if (pos[i] < rect->min[i]) {
			result += SQ(rect->min[i] - pos[i]);
		} else if (pos[i] > rect->max[i]) {
			result += SQ(rect->max[i] - pos[i]);
		}
	}

	return result;
}




static kdres* kd_hp_create(int size)
{
    kdres* heap  = (kdres *)malloc(sizeof(kdres));
    heap->data = (res_node*)malloc(size*sizeof(res_node));
    heap->size = 0;


    return heap;
}

static void kd_hp_swap(res_node *node1, res_node *node2)
{
    res_node tmp;
    tmp = *node1;
    *node1 = *node2;
    *node2 = tmp;

}

/*

 Heapify function is used to make sure that the heap property is never violated .
 In case of deletion of a node, or creating a max heap form aan array, heap propety
may be violated. In such cases, Heaprity function can be called to make sure that heap
propety is never violated.

*/
static void kd_heapify(kdres *heap, int i)
{
    int largest;
    if(i < heap->size)
    {
        largest = (HP_LCHILD(i) < heap->size && heap->data[HP_LCHILD(i)].dist_sq > heap->data[i].dist_sq)?HP_LCHILD(i):i;
        if(HP_RCHILD(i) < heap->size && heap->data[HP_RCHILD(i)].dist_sq > heap->data[largest].dist_sq)
        {
            largest = HP_RCHILD(i);
        }
        if(largest != i)
        {
            kd_hp_swap(&(heap->data[i]), &(heap->data[largest]));
            kd_heapify(heap, largest);
        }

    }

}
static int kd_hp_insert(kdres *heap, res_node node)
{
    int i = (heap->size)++;
    while(i && node.dist_sq > heap->data[HP_PARENT(i)].dist_sq)
    {
        heap->data[i] = heap->data[HP_PARENT(i)];
        i = HP_PARENT(i);

    }
    heap->data[i] = node;
    return 1;

}

static int kd_hp_modify(kdres *heap, res_node node)
{
    int i;

    heap->data[0] = node;
    kd_heapify(heap, 0);
    return 1;

    return -1;
}


static void kd_hp_sort(kdres *heap)
{
    int size = heap->size;
    int i = size - 1;
    while(i > 0)
    {
        kd_hp_swap(&(heap->data[0]), &(heap->data[i]));
        --heap->size;
        kd_heapify(heap, 0);
        --i;
    }
    heap->size = size;
}
//以遞歸結構回溯查找
static int find_nearest_n(struct kdnode *node, const double *pos, kdres *result, struct kdhyperrect* rect, int num)
{
    int i;
    double dummy, dist_sq;
    struct kdnode *nearer_subtree, *farther_subtree;
    double *nearer_hyperrect_coord, *farther_hyperrect_coord;
    res_node hp_node;

    if(!node)
    {
        return 0;
    }

    int dir = node->dir;
    /* Check the distance of the point at the current node, compare it
     * with our best so far */
    dist_sq = 0;
    for(i=0; i < rect->dim; i++)
    {
        dist_sq += SQ(node->pos[i] - pos[i]);
    }

    hp_node.item = node;
    hp_node.dist_sq = dist_sq;
    //如果堆的大小大于num,也就是大于总的要找的节点数

    if(result->size == num && result->data[0].dist_sq > dist_sq)
    {
        /* get furthest element */
        //得到最远的节点
        kd_hp_modify(result, hp_node);
    }
    else if(result->size < num)
    {
//        如果堆的大小小于num,直接将此节点插入堆中
        if(kd_hp_insert(result, hp_node) == -1)
        {
            return 0;
        }
    }
    //在二叉树中,决定向左走还是向右走
    /* Decide whether to go left or right in the tree */
    dummy = pos[dir] - node->pos[dir];
    if (dummy <= 0)
    {
        nearer_subtree = node->left;
        farther_subtree = node->right;
        nearer_hyperrect_coord = rect->max + dir;
        farther_hyperrect_coord = rect->min + dir;
    }
    else
    {
        nearer_subtree = node->right;
        farther_subtree = node->left;
        nearer_hyperrect_coord = rect->min + dir;
        farther_hyperrect_coord = rect->max + dir;
    }

    if (nearer_subtree)
    {
        /* Slice the hyperrect to get the hyperrect of the nearer subtree */
        dummy = *nearer_hyperrect_coord;
        *nearer_hyperrect_coord = node->pos[dir];
        /* Recurse down into nearer subtree */
        find_nearest_n(nearer_subtree, pos, result, rect,num);
        /* Undo the slice */
        *nearer_hyperrect_coord = dummy;
    }


    if (farther_subtree)
    {
        /* Get the hyperrect of the farther subtree */
        dummy = *farther_hyperrect_coord;
        *farther_hyperrect_coord = node->pos[dir];
        /* Check if we have to recurse down by calculating the closest
         * point of the hyperrect and see if it's closer than our
         * minimum distance in result_dist_sq. */
       // if (hyperrect_dist_sq_dir(rect, pos,*farther_hyperrect_coord)  < rheap_get_max(*result)->dist_sq) {
        if (hyperrect_dist_sq(rect, pos) < result->data[0].dist_sq)
        {
            /* Recurse down into farther subtree */
            find_nearest_n(farther_subtree, pos, result, rect,num);
        }
        /* Undo the slice on the hyperrect */
        *farther_hyperrect_coord = dummy;
    }
    return 1;
}
//@pos target position
kdres* kd_nearest_n(struct kdtree *kd, const double *pos,int num)
{
    struct kdhyperrect *rect;
    kdres *rset;

    if (!kd) return 0;
    if (!kd->rect) return 0;

    /* Allocate result set */
    if(!(rset = (kdres*)malloc(sizeof(kdres)))) {
        return 0;
    }


//复制边界超矩形，我们将处理副本
    /* Duplicate the bounding hyperrectangle, we will work on the copy */
    if (!(rect = hyperrect_duplicate(kd->rect))) {
        kd_res_free(rset);
        return 0;
    }
    rset->data = (res_node*)malloc(num*sizeof(res_node));
    rset->size=0;
//    rset->max_size = num;

//我们的第一个猜测是根节点
    /* Our first guesstimate is the root node */
    /* Search for the nearest neighbour recursively */
    find_nearest_n(kd->root, pos, rset, rect,num);

//    kd_hp_sort(rset);
    /* Free the copy of the hyperrect */
    hyperrect_free(rect);

    return rset;

}

static int find_nearest_range(struct kdnode *node, const double *pos, double range, kdres* result, int dim)
{

    double dist_sq, dx;
    int i, ret, added_res = 0;
    res_node hp_node;
    if(!node)
    {
        return 0;
    }

    int dir = node->dir;
    /* if the photon is close enough, add it to the result heap */
    dist_sq = 0;
    for(i=0; i<dim; i++) {
        dist_sq += SQ(node->pos[i] - pos[i]);
    }

    if(dist_sq <= range*range)
    {
        hp_node.item = node;
        hp_node.dist_sq = dist_sq;
        if(result->size)
        {
            result->data = (res_node*)realloc(result->data, (result->size + 1) * sizeof(res_node));

            kd_hp_insert(result, hp_node);
        }
        else
        {
            result->data = (res_node*)malloc(sizeof (res_node));
            kd_hp_insert(result, hp_node);

        }

        added_res = 1;
    }
//	/* find signed distance from the splitting plane */
    dx = pos[dir] - node->pos[dir];

    ret = find_nearest_range(dx <= 0.0 ? node->left : node->right, pos, range, result, dim);
    if(ret >= 0 && fabs(dx) < range)
    {
        added_res += ret;
        ret = find_nearest_range(dx <= 0.0 ? node->right : node->left, pos, range, result, dim);
    }

}


kdres *kd_nearest_range(struct kdtree *kd, const double *pos, double range)
{
    kdres *rset;
    if(!(rset = (kdres *)malloc(sizeof *rset))) {
        return 0;
    }

    rset->size=0;
    rset->data = NULL;

    find_nearest_range(kd->root, pos, range, rset, kd->dim);
    kd_hp_sort(rset);
    return rset;
}
kdres *kd_nearest_rangef(struct kdtree *kd, const float *pos, float range)
{
	static double sbuf[16];
	double *bptr, *buf = 0;
	int dim = kd->dim;
	kdres *res;

	if(dim > 16) {
#ifndef NO_ALLOCA
		if(dim <= 256)
			bptr = buf = (double*)alloca(dim * sizeof *bptr);
		else
#endif
			if(!(bptr = buf = (double*)malloc(dim * sizeof *bptr))) {
				return 0;
			}
	} else {
		bptr = buf = sbuf;
	}

	while(dim-- > 0) {
		*bptr++ = *pos++;
	}

	res = kd_nearest_range(kd, buf, range);

#ifndef NO_ALLOCA
	if(kd->dim > 256)
#else
	if(kd->dim > 16)
#endif
		free(buf);
	return res;
}

kdres *kd_nearest_range3(struct kdtree *tree, double x, double y, double z, double range)
{
	double buf[3];
	buf[0] = x;
	buf[1] = y;
	buf[2] = z;
	return kd_nearest_range(tree, buf, range);
}
kdres *kd_nearest_range3f(struct kdtree *tree, float x, float y, float z, float range)
{
	double buf[3];
	buf[0] = x;
	buf[1] = y;
	buf[2] = z;
	return kd_nearest_range(tree, buf, range);
}
