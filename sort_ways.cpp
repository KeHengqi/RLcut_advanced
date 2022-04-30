#include "sort_ways.hpp"
#include <stdlib.h>

// a for array, n for number of array
// Ascendent array
// Checked
void Bubble_sort(double* a, int n)
{
	int i;
	for (int i = 0; i < n - 1; i++)
	{
		int sign = 0;
		for (int j = 0; j < n - i - 1; j++)
		{
			if (a[j] > a[j + 1])
			{
				double temp = a[j];
				a[j] = a[j + 1];
				a[j + 1] = temp;
				sign = 1;
			}
		}
		if (sign == 0)
		{
			break;
		}
	}
}
// Checked
void Selection_sort(double* a, int n)
{
	for (int i = 0; i < n; i++)
	{
		int min_index = i;
		for (int j = i; j < n; j++)
		{
			if (a[min_index] > a[j])
			{
				min_index = j;
			}
		}
		if (min_index != i)
		{
			double temp = a[min_index];
			a[min_index] = a[i];
			a[i] = temp;
		}
	}
}

// Checked
void Insert_sort(double* a, int n)
{
	for (int i = 1; i < n; i++)
	{
		for (int j = i; j > 0; j--)
		{
			if (a[j] < a[j - 1])
			{
				double temp = a[j];
				a[j] = a[j - 1];
				a[j - 1] = temp;
			}
			else
			{
				break;
			}
		}
	}
}

// Checked
void quick_sort_real(double* start, double* end, double standard)
{
	double* work_p_left = start;
	double* work_p_right = end;
	while (work_p_left < work_p_right)
	{
		while (*work_p_right > standard)
		{
			work_p_right--;
		}
		double temp = *work_p_right;
		*work_p_right = *work_p_left;
		*work_p_left = temp;
		while (*work_p_left < standard)
		{
			work_p_left++;
		}
		temp = *work_p_right;
		*work_p_right = *work_p_left;
		*work_p_left = temp;
	}
	*work_p_left = standard;
	if (work_p_left != start)
	{
		quick_sort_real(start, work_p_left - 1, *start);
	}
	if (work_p_right != end)
	{
		quick_sort_real(work_p_right + 1, end, *(work_p_right + 1));
	}
}
void Quick_sort(double* a, int n){
	quick_sort_real(a, a + n - 1, a[0]);
};

void merge_sort_real(double* array1, double* array2, int a1_size, int a2_size)
{
	if (a1_size == 2)
	{
		if (*array1 > *(array1 + 1))
		{
			double temp = *array1;
			*array1 = *(array1 + 1);
			*(array1 + 1) = temp;
		}
	}
	if (a2_size == 2)
	{
		if (*array2 > *(array2 + 1))
		{
			double temp = *array2;
			*array2 = *(array2 + 1);
			*(array2 + 1) = temp;
		}
	}

	if (a1_size > 2)
	{
		int half = a1_size / 2;
		merge_sort_real(array1, array1 + half, half, a1_size - half);
	}
	if (a2_size > 2)
	{
		int half = a2_size / 2;
		merge_sort_real(array2, array2 + half, half, a2_size - half);
	}

	double* temp_array = (double*)malloc(sizeof(double) * (a1_size + a2_size));
	double* temp_work = temp_array;
	int a1_work, a2_work;
	//double* a1_work, * a2_work;
	a1_work = 0;
	a2_work = 0;
	while (a1_work < a1_size && a2_work < a2_size)
	{
		if (*(array1 + a1_work) < *(array2 + a2_work))
		{
			*temp_work = *(array1 + a1_work);
			//temp_work++;
			a1_work++;
		}
		else
		{
			*temp_work = *(array2 + a2_work);
			a2_work++;
		}
		temp_work++;
	}
	while (a1_work < a1_size)
	{
		*temp_work = *(array1 + a1_work);
		a1_work++;
		temp_work++;
	}
	while (a2_work < a2_size)
	{
		*temp_work = *(array2 + a2_work);
		a2_work++;
		temp_work++;
	}

	for (a1_work = 0; a1_work < a1_size + a2_size; a1_work++)
	{
		*(array1 + a1_work) = temp_array[a1_work];
	}
	free(temp_array);

}
void Merge_sort(double* a, int n)
{
	double* a1 = a;
	int half = n / 2;
	double* a2 = a + half;
	merge_sort_real(a1, a2, half, n - half);
}

// Checked
void max_heaptify(double* tree_head, int heap_head_index, int size_of_tree)
{
	int min_non_tree_index = size_of_tree / 2 - 1;
	int left_child = 2 * heap_head_index + 1;
	int right_child = 2 * heap_head_index + 2;
	int max_child = -1;
	if (heap_head_index > min_non_tree_index)
	{
		return;
	}
	if (right_child >= size_of_tree)
	{
		max_child = left_child;
	}
	else
	{
		if (tree_head[left_child] > tree_head[right_child])
		{
			max_child = left_child;
		}
		else
		{
			max_child = right_child;
		}
	}
	if (tree_head[max_child] > tree_head[heap_head_index])
	{
		double temp = tree_head[max_child];
		tree_head[max_child] = tree_head[heap_head_index];
		tree_head[heap_head_index] = temp;
		max_heaptify(tree_head, max_child, size_of_tree);
	}
}
void max_heaptify_int(int* tree_head, int heap_head_index, int size_of_tree)
{
    int min_non_tree_index = size_of_tree / 2 - 1;
	int left_child = 2 * heap_head_index + 1;
	int right_child = 2 * heap_head_index + 2;
	int max_child = -1;
	if (heap_head_index > min_non_tree_index)
	{
		return;
	}
	if (right_child >= size_of_tree)
	{
		max_child = left_child;
	}
	else
	{
		if (tree_head[left_child] > tree_head[right_child])
		{
			max_child = left_child;
		}
		else
		{
			max_child = right_child;
		}
	}
	if (tree_head[max_child] > tree_head[heap_head_index])
	{
		int temp = tree_head[max_child];
		tree_head[max_child] = tree_head[heap_head_index];
		tree_head[heap_head_index] = temp;
		max_heaptify_int(tree_head, max_child, size_of_tree);
	}
}
void max_heaptify_long(long* tree_head, int heap_head_index, int size_of_tree)
{
    int min_non_tree_index = size_of_tree / 2 - 1;
	int left_child = 2 * heap_head_index + 1;
	int right_child = 2 * heap_head_index + 2;
	int max_child = -1;
	if (heap_head_index > min_non_tree_index)
	{
		return;
	}
	if (right_child >= size_of_tree)
	{
		max_child = left_child;
	}
	else
	{
		if (tree_head[left_child] > tree_head[right_child])
		{
			max_child = left_child;
		}
		else
		{
			max_child = right_child;
		}
	}
	if (tree_head[max_child] > tree_head[heap_head_index])
	{
		long temp = tree_head[max_child];
		tree_head[max_child] = tree_head[heap_head_index];
		tree_head[heap_head_index] = temp;
		max_heaptify_long(tree_head, max_child, size_of_tree);
	}   
}
void max_heaptify_ll(long long* tree_head, int heap_head_index, int size_of_tree)
{
    int min_non_tree_index = size_of_tree / 2 - 1;
	int left_child = 2 * heap_head_index + 1;
	int right_child = 2 * heap_head_index + 2;
	int max_child = -1;
	if (heap_head_index > min_non_tree_index)
	{
		return;
	}
	if (right_child >= size_of_tree)
	{
		max_child = left_child;
	}
	else
	{
		if (tree_head[left_child] > tree_head[right_child])
		{
			max_child = left_child;
		}
		else
		{
			max_child = right_child;
		}
	}
	if (tree_head[max_child] > tree_head[heap_head_index])
	{
		long long temp = tree_head[max_child];
		tree_head[max_child] = tree_head[heap_head_index];
		tree_head[heap_head_index] = temp;
		max_heaptify_ll(tree_head, max_child, size_of_tree);
	}
}
void min_heaptify(double* tree_head, int heap_head_index, int size_of_tree)
{
	int min_non_tree_index = size_of_tree / 2 - 1;
	int left_child = 2 * heap_head_index + 1;
	int right_child = 2 * heap_head_index + 2;
	int min_child = -1;
	if (heap_head_index > min_non_tree_index)
	{
		return;
	}
	if (right_child >= size_of_tree)
	{
		min_child = left_child;
	}
	else
	{
		if (tree_head[left_child] < tree_head[right_child])
		{
			min_child = left_child;
		}
		else
		{
			min_child = right_child;
		}
	}
	if (tree_head[min_child] < tree_head[heap_head_index])
	{
		double temp = tree_head[min_child];
		tree_head[min_child] = tree_head[heap_head_index];
		tree_head[heap_head_index] = temp;
		min_heaptify(tree_head, min_child, size_of_tree);
	}
}
void Hash_sort(double* a, int n)
{
	//n++;
	int min_non_tree_index = n / 2 - 1;
	for (int i = min_non_tree_index; i >= 0; i--)
	{
		max_heaptify(a, i, n);
	}
	for (int j = n - 1; j >= 0; j--)
	{
		double temp = a[0];
		a[0] = a[j];
		a[j] = temp;
		max_heaptify(a, 0, j);
	}
};

