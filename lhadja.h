#define __LHADJA_H

#ifndef __LHADJA_H
class ART
{
public:
    ART():value = -1;
    int getValue() const { return this->value; } // définition au sein de la classe
    void setValue(int valuer);
    void print() const;
private:
    int value;
};

#endif