#include "utils.h"

//    template <typename ValueType> 
    void rotationplan(ValueType& dx, ValueType& dy, ValueType& cs, ValueType& sn)
    {
      ValueType temp = cs * dx + sn *dy;
      dy = -sn*dx+cs*dy;
      dx = temp;
    }

//    template <typename ValueType>
    void genererrotaionplan(ValueType& dx, ValueType& dy, ValueType& cs, ValueType& sn)
    {
      if(dy == ValueType(0.0)){
			cs = 1.0;
			sn = 0.0;
      }else if (abs(dy) > abs(dx)) {
			ValueType tmp = dx / dy;
			sn = ValueType(1.0) / sqrt(ValueType(1.0) + tmp*tmp);
			cs = tmp*sn;            
      }else {
			ValueType tmp = dy / dx;
			cs = ValueType(1.0) / sqrt(ValueType(1.0) + tmp*tmp);
			sn = tmp*cs;
      }
    }

//    template <class LinearOperator,typename ValueType> 
    void applyrotationplan(LinearOperator& H, ValueType& cs, ValueType& sn, ValueType& s, int i)
    {
      for (int k = 0; k < i; k++){
			rotationplan(H(k,i), H(k+1,i), cs[k], sn[k]);
      }
      GeneratePlaneRotation(H(i,i), H(i+1,i), cs[i], sn[i]);
      rotationplan(H(i,i), H(i+1,i), cs[i], sn[i]);
      rotationplan(s[i], s[i+1], cs[i], sn[i]);
    }
