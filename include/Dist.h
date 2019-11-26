#ifndef DIST_H
#define DIST_H

#include <algorithm>
#include <vector>
#include "EuclideanDistance.h"

using namespace std;

FD_NS_START
class Dist
{
    public:
        Dist(){};
        virtual ~Dist(){};

        template <typename ValueType,JInt nDimension>
        ValueType calcDistance(const MeasurementVector<ValueType, nDimension>& v1, const MeasurementVector<ValueType, nDimension>& v2) const
        {
            ValueType val;
            double exp = 0.5;
            double mi = min(v1[0], v2[0])+exp;
            double ma = max(v1[0], v2[0])+exp;
            val = (ma/mi-1)*max(v1[1], v2[1]);
            return val;
        }

        template <typename ValueType>
        ValueType calcDistance(const std::vector<ValueType>& v1, const std::vector<ValueType>& v2) const
        {
            ValueType val;
            double exp = 0.5;
            double mi = min(v1[0], v2[0])+exp;
            double ma = max(v1[0], v2[0])+exp;
            val = (ma/mi-1)*max(v1[1], v2[1]);
            return val;
        }

    protected:

    private:
};
FD_NS_END

#endif // DIST_H
