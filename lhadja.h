#define __LHADJA_H

#ifndef __LHADJA_H
class A
{
public:
    A():value = -1;
    int getValue() const { return this->value; } // définition au sein de la classe
    void setValue(int value);
    void print() const;
private:
    int value;
};

#endif