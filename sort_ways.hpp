#pragma once
#ifndef SORT_WAYS_H
#define SORT_WAYS_H
#endif // !SORT_WAYS_H

typedef void(*Sort_method)(double*, int);

void Selection_sort(double* a, int n);
void Bubble_sort(double* a, int n);
void Bubble_sort_advanced(double* a, int n);
void Insert_sort(double* a, int n);
void Quick_sort(double* a, int n);
void max_heaptify(double* tree_head, int heap_head_index, int size_of_tree);
void max_heaptify_int(int* tree_head, int heap_head_index, int size_of_tree);
void max_heaptify_long(long* tree_head, int heap_head_index, int size_of_tree);
void max_heaptify_ll(long long* tree_head, int heap_head_index, int size_of_tree);
void min_heaptify(double* tree_head, int heap_head_index, int size_of_tree);
void Hash_sort(double* a, int n);
void Merge_sort(double* a, int n);
int Test_sort(Sort_method test_function, unsigned int seed);
