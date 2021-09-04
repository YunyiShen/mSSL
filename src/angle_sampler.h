#ifndef ANGLE_SAMPLER_H
#define ANGLE_SAMPLER_H
#include "active_intersections.h"

class AngleSampler
{
public:
    ActiveIntersections active_intersections;
    double rotation_angle;
    arma::mat rotated_slices;
    AngleSampler() = default;
    AngleSampler(ActiveIntersections active_intersect)
    {
        active_intersections = active_intersect;
        rotationAngle temp = active_intersect.rotated_intersections();
        rotation_angle = temp.rotation_angle;
        arma::vec slices = temp.shifted_angle;
        rotated_slices = slices;
        rotated_slices.reshape(slices.n_elem / 2, 2);
    }

    arma::vec _get_slices_cumulative_length();
    double draw_angle();
};

arma::vec AngleSampler::_get_slices_cumulative_length()
{
    arma::vec cum_len(rotated_slices.n_rows + 1, fill::zeros);
    for (int i = 0; i < rotated_slices.n_rows; i++)
    {
        cum_len(i + 1) = cum_len(i) + rotated_slices(i, 1) - rotated_slices(i, 0);
    }
    return (cum_len);
}

double AngleSampler::draw_angle()
{
    arma::vec cum_len = this->_get_slices_cumulative_length();
    double l = cum_len(cum_len.n_elem - 1);
    double sample = l * arma::randu<double>();
    int idx = 1;
    for (idx = 1; idx < cum_len.n_elem; idx++)
    {
        if (cum_len(idx - 1) < sample && cum_len(idx) >= sample)
        {
            break;
        }
    }
    idx--;
    return (rotated_slices(idx, 0) + sample - cum_len(idx) + rotation_angle);
}

#endif