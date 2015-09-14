#include "lhadja.h"

inline void ART::setValue(int valuer) // définition en dehors de la classe (inline)
{
    value = valuer;
}

void ART::print() const // définition en dehors de la classe (non inline)
{
    std::cout << "Value=" << value << std::endl;
}