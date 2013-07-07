/*************************************************************************************
 *
 * MATLAB (R) is a trademark of The Mathworks (R) Corporation
 *
 * Filename:    mxGetPropertyPtr.c
 * Programmer:  James Tursa
 * Version:     1.00
 * Date:        March 07, 2011
 * Copyright:   (c) 2011 by James Tursa, All Rights Reserved
 *
 *  This code uses the BSD License:
 *
 *  Redistribution and use in source and binary forms, with or without 
 *  modification, are permitted provided that the following conditions are 
 *  met:
 *
 *     * Redistributions of source code must retain the above copyright 
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright 
 *       notice, this list of conditions and the following disclaimer in 
 *       the documentation and/or other materials provided with the distribution
 *      
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *-------------------------------------------------------------------------
 *
 *  This file is intended to be included in user source code to provide the
 *  capability to get the pointer to a property inside a C mex routine. The
 *  signature is as follows:
 *
 *  #include "matrix.h"
 *  mxArray *mxGetPropertyPtr(const mxArray *pa, mwIndex index, const char *propname);
 *
 *  If successful, the returned value will be a pointer to the named property
 *  of the object. Note that this differs from the behavior of the API routine
 *  mxGetProperty:
 *
 *    mxGetPropertyPtr: Returns a pointer to the actual property itself.
 *    mxGetProperty:    Returns a pointer to a *copy* of the property.
 *
 *  If unsuccessful, the mxGetPropertyPtr will return NULL.
 *
 *  Limitations of mxGetPropertyPtr:
 *  - Works on classdef objects, which require MATLAB R2008a or later.
 *  - Generally, the property that this pointer points to should be treated
 *    as read-only unless you REALLY, REALLY, KNOW WHAT YOU ARE DOING!
 *  - Works only in mex routines, not in Engine applications.
 *  - Requires a helper C-mex routine, GetPropertyPtrx.
 *  - Requires a helper m-file routine, GetPropertyPtr.
 *  - mxGetPropertyPtr uses unofficial behind-the-scenes hacking of the
 *    mxArray structure itself to recover the pointer. There is no guarantee
 *    that this will work on all versions of MATLAB. If it doesn't work for
 *    you, please contact the author with details of your problem.
 *
 *  Although the intended use of mxGetPropertyPtr is for the newer classdef
 *  class variables (for which TMW did not supply an equivalent routine), in
 *  fact mxGetPropertyPtr will also work with structures and old-style user-
 *  define classes as well, treating propname as a fieldname.
 *
 *  Using mxGetPropertyPtr is much more efficient for large size properties
 *  than using mxGetProperty because no large intermediate copy of the
 *  property is created. This will result in speed increases as well as
 *  avoiding potential out-of-memory conditions that can result from using
 *  mxGetProperty. Also, on some versions of MATLAB the mxGetProperty function
 *  can return an invalid pointer in out-of-memory cases and using this
 *  pointer will crash MATLAB. For small size properties, the extra overhead
 *  involved with mxGetPropertyPtr (calling mexCallMATLAB with calls the
 *  m-file GetPropertyPtr which calls another mex routine GetPropertyPtrx)
 *  outweighs its benefits and one would be better off calling the API
 *  function mxGetProperty instead. But if you don't know the size of the
 *  property in advance and the size could be large, it would be advisable
 *  to use mxGetPropertyPtr.
 *
 * Change Log:
 * 2011/Mar/07 --> 1.00, Initial Release
 *
 ****************************************************************************/

#ifndef  mxGetPropertyPtr_c  /* Include guard */
#define  mxGetPeopertyPtr_c

/* Includes ----------------------------------------------------------- */

#include "mex.h"

/* Macros ------------------------------------------------------------- */

#ifndef  INT64_TYPE
#define  INT64_TYPE long long
#endif

#ifndef  MWSIZE_MAX
#define  MWSIZE_MAX
#define  mwIndex        int
#define  mwSignedIndex  int
#define  mwSize         int
#endif

/* -------------------------------------------------------------------- *
 * mexCallMATLABWithTrap first introduced in R2008b.                    *
 * For earlier versions use this intermediate interface instead.        *
 * The macro mxSetLogical is set in matrix.h in R2008a and earlier, but *
 * not in R2008b and later. Use that to determine need for trap code.   *
 * -------------------------------------------------------------------- */

#ifdef  mxSetLogical
mxArray *mexCallMATLABWithTrap(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], const char *functionName)
{
	mexSetTrapFlag(1);
	if( mexCallMATLAB(nlhs,plhs,nrhs,prhs,functionName) ) {
		return mxCreateDoubleScalar(0.0); /* Return something (anything) that is non-NULL to indicate error */
	} else {
		return NULL;
	}
}
#endif

/* -------------------------------------------------------------------- */

mxArray *mxGetPropertyPtr(const mxArray *pa, mwIndex ix, char *propname)
{
	mxArray *mx, *PropertyPtr = NULL;
	mxArray *lhs[1], *rhs[3];
	union {INT64_TYPE theinteger; mxArray *thepointer;} uip; /* For non-compliant code use below */

	if( pa && propname && ix >= 0 ) {
        if( mxIsStruct(pa) ) {  /* If a struct, just use mxGetField directly */
            PropertyPtr = mxGetField(pa,ix,propname);
        } else {  /* Otherwise, assume a classdef object and go get the pointer */
		    rhs[0] = (mxArray *) pa;
		    rhs[1] = mxCreateDoubleScalar(ix+1);
		    rhs[2] = mxCreateString(propname);
			mx = mexCallMATLABWithTrap(1,lhs,3,rhs,"GetPropertyPtr");
		    if( !mx && lhs[0] ) {
		    	uip.theinteger = *((INT64_TYPE *)mxGetData(lhs[0]));
		    	PropertyPtr = uip.thepointer; /* Need to use non-compliant union code to avoid compiler bug */
		    	mxDestroyArray(lhs[0]);
		    }
			if( mx ) mxDestroyArray(mx);
		    mxDestroyArray(rhs[2]);
		    mxDestroyArray(rhs[1]);
        }
	}
	return PropertyPtr;
}

#endif  /* mxGetPropertyPtr_c */
