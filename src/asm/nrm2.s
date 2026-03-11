## nrm2(n, x) -- 2-norm
## returns: |x|_2


## float (SSE single)
.global blas_nrm2f
blas_nrm2f:


## TODO

## NOTE: sqrtps
ret


## double (SSE double)
.global blas_nrm2
.type blas_nrm2, @function
blas_nrm2:

## TODO

## NOTE: sqrtpd
ret


## long double (x87 extended floating point)
.global blas_nrm2l
.type blas_nrm2l, @function
blas_nrm2l:

## TODO

## NOTE: fsqrt 
ret

## vim: set ft=gas commentstring=##\ %s vartabstop=2,10,25,4 noexpandtab:
