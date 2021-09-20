#ifndef QUIC_UNSHORT_PAIR_H
#define QUIC_UNSHORT_PAIR_H

#ifndef EPS
#define EPS (double(2.22E-16))
#endif

namespace quic{

typedef struct
{
	int i;
	int j;
} ushort_pair_t;

}
#endif