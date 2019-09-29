#pragma once
#include <fixed_point.h>
#include <privacyconf.h>
#include <iostream>
#include <string>

template <typename T>
void printVector(std::string s, T* v, int size) {
  printf("%s\n", s.c_str());
  for (int i = 0; i < size; i++) {
    if (printints)
      printf("%d", *(int*)&v[i]);
    else {
      // if (typeid(T) == typeid(float) || typeid(T) == typeid(double)) {
      //     printf("%.8lf ", v[i]);
      // } else if (typeid(T) == typeid(fixed_point<int64_t,32>)) {
      //     //std::cout.precision(8);
      //     std::cout << v[i] << " ";
      // } else if (typeid(T) == typeid(fixed_point<int32_t,16>)) {
      //     //std::cout.precision(8);
      //     std::cout << v[i] << " ";
      // } else {
      //     std::cout << v[i] << " ";
      // }
      std::cout << v[i] << " ";
    }
  }
  printf("\n\n");
}

template <typename T>
void printMatrix(std::string s, T** a, int rows, int cols) {
  printf("%s\n", s.c_str());
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      if (printints)
        printf("%d", *(int*)&a[i][j]);
      else {
        // if (typeid(T) == typeid(float) || typeid(T) == typeid(double)) {
        //     printf("%.8lf ", a[i][j]);
        // } else if (typeid(T) == typeid(fixed_point<int64_t,32>)) {
        //     //std::cout.precision(8);
        //     std::cout << a[i][j] << " ";
        // } else if (typeid(T) == typeid(fixed_point<int32_t,16>)) {
        //     //std::cout.precision(8);
        //     std::cout << a[i][j] << " ";
        // } else {
        //     std::cout << a[i][j] << " ";
        // }
        std::cout << a[i][j] << " ";
      }
    }
    printf("\n");
  }
  printf("\n");
}
