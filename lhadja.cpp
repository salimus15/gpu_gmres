#include "lhadja.h"

inline void A::setValue(int value) // définition en dehors de la classe (inline)
{
    this->value = value;
}

void A::print() const // définition en dehors de la classe (non inline)
{
    std::cout << "Value=" << this->value << std::endl;
}