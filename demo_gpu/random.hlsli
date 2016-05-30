#ifndef __RANDOM__
#define __RANDOM__

void psdes(inout uint lword, inout uint irword)
{
	const uint NITER=4;
	const uint c1[NITER]={
		0xbaa96887L, 0x1e17d32cL, 0x03bcdc3cL, 0x0f33d1b2L};
	const uint c2[NITER]={
		0x4b0f3b58L, 0xe874f0c3L, 0x6955c5a6L, 0x55a7ca46L};
	uint i,ia,ib,iswap,itmph=0,itmpl=0;

	for (i=0;i<NITER;i++) {
		ia=(iswap=irword) ^ c1[i];
		itmpl = ia & 0xffff;
		itmph = ia >> 16;
		ib=itmpl*itmpl+ ~(itmph*itmph);
		irword=lword ^ (((ia = (ib >> 16) |
			((ib & 0xffff) << 16)) ^ c2[i])+itmpl*itmph);
		lword=iswap;
	}
}

// Generate random float in [0, 1)
float rnd(inout int2 idum)
{
	const uint jflone = 0x3f800000;
	const uint jflmsk = 0x007fffff;
	uint irword,itemp,lword;
	
	if (idum.x < 0) {
		idum.y = -idum.x;
		idum.x=1;
	}
	irword=idum.x;
	lword=idum.y;
	psdes(lword,irword);
	itemp=jflone | (jflmsk & irword);
	++idum.x;
	return asfloat(itemp)-1.0f;
}

#endif