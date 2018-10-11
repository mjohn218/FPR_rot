#include "vector_rot_calls.h"

void rotationEuler(double tx, double ty, double tz, double *M) {
	/*Using yaw, pitch, roll rotation matrices.
	 Extrinsic rotations around x, then y, then z. 
	 */
	/*Rx=[1 0 0; 0 costhet -sinthet ; 0 sinthet costhet]
	 Ry=[costhet 0 sinthet; 0 1 0 ; -sinthet 0 costhet]
	 Rz=[costhet -sinthet 0 ; sinthet costhet 0 ; 0 0 1]*/

	/*Ryx=[ cosy sinxsiny sinycosx ;
	 0        cosx -sinx ; 
	 -siny cosysinx cosycosx]*/

	/*Rzyx=[  coszcosy, coszsinxsiny -sinzcosx, coszsinycosx+sinzsinx;
	 sinzcosy, sinzsinxsiny+coszcosx, sinzsinycosx-coszsinx;
	 -siny, cosysinx, cosycosx]
	 */

	double sx = sin(tx);
	double cx = cos(tx);
	double sy = sin(ty);
	double cy = cos(ty);
	double sz = sin(tz);
	double cz = cos(tz);
	M[0] = cz * cy;
	M[1] = cz * sx * sy - sz * cx;
	M[2] = cz * sy * cx + sz * sx;
	M[3] = sz * cy;
	M[4] = sz * sx * sy + cz * cx;
	M[5] = sz * sy * cx - cz * sx;
	M[6] = -sy;
	M[7] = cy * sx;
	M[8] = cy * cx;

}
