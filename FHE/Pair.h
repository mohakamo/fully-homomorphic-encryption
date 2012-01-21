/* Pair.h
 * Created by valerini at 6/1/2012
 */
#ifndef _PAIR_H_
#define _PAIR_H_

template<typename FirstType, typename SecondType>
  class Pair {
 public:
  FirstType first;
  SecondType second;
  /* public:
  FirstType Get_First() {
    return first;
  }

  SecondType Get_Second() {
    return second;
    }*/

  Pair(FirstType first, SecondType second) {
    this->first = first;
    this-> second = second;
  }
};

#endif /* _PAIR_H_ */
