/*
 queue.cpp
 Created: 2/14/2012
 Author: T. Bui

*/

#include "queue.h"
#include <cassert>
#include <iostream>
using namespace std;

Queue::Queue(int c)
{
  assert(c > 0);
  qCapacity = c;
  qArray = new int[qCapacity];
  front = -1;
  back = 0;
  qLength = 0;
}//Queue

Queue::Queue(const Queue& source)
{
  qCapacity = source.qCapacity;
  qArray = new int[qCapacity];
  front = source.front;
  back = source.back;
  qLength = source.qLength;
  for (int i = 0; i < qCapacity; i++)
    qArray[i] = source.qArray[i];
}//copy constructor

Queue::~Queue()
{
  delete [] qArray;
}

void Queue::enqueue(int key)
{
  assert(qLength < qCapacity);
  qArray[back] = key;
  cout << "enqueue " << key << " to cell " << back << endl;
  if (qLength == 0)  //special case when queue was empty
    front = back;    //need to initialize front properly
  back = (back + 1) % qCapacity;
  qLength++;
}//enqueue

int Queue::dequeue()
{
  assert(!isEmpty());
  int x = qArray[front];
  cout << "dequeue " << x << " from cell " << front << endl;
  if (qLength == 1)  //queue is empty after dequeue
    front = -1;      //update front appropriately
  else 
    front = (front + 1) % qCapacity;
  qLength--;
  return x;
}//dequeue
  
void Queue::clear()
{
  front = -1;
  back = 0;
  qLength = 0;
}//clear
