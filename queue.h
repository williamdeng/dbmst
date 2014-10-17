/*
  queue.h
  Date: 2/14/2012
  Author: T. Bui
*/


#ifndef QUEUE_H
#define QUEUE_H
#include <cassert>

class Queue{
private:
  int* qArray;
  int qLength;   //current number of elements in the queue
  int qCapacity;  // capacity of the queue
  int front;  //points to the first element in the queue
  int back;   //points to the next available position in the queue


public:
  Queue(int size);
  Queue(const Queue& source);
  ~Queue();
  void enqueue(int key);
  int dequeue();
  int getLength() const{return qLength;}
  int getCapacity() const{return qCapacity;}
  bool isEmpty() const{return (qLength == 0);}
  bool isFull() const{return (qLength == qCapacity);}
  void clear();
};//Queue

#endif
