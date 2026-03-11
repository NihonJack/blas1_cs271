## axpyf(n, a, x, y) -- update vector
## returns: y = a * x + y


## NOTE: Pay close attention to arguments
## 	n (INT) -> rdi
##	a (SSE) -> xmm0 / ST(0)
##	x (INT) -> rsi
##	y (INT) -> rdx


## float (SSE single)
.global blas_axpyf
.type blas_axpyf, @function
blas_axpyf:
## NOTE: xmm0 contains alpha

## NOTE: Use the appropriate fused multiply-add instruction
##	vfmadd132ps/vfmadd213ps/vfmadd231ps

## TODO

ret


## double (SSE double)
.global blas_axpy
.type blas_axpy, @function
blas_axpy:
## NOTE: xmm0 contains alpha

## NOTE: Use the appropriate fused multiply-add instruction
##	vfmadd132pd/vfmadd213pd/vfmadd231pd

## TODO

ret


## long double (x87 extended floating point)
.global blas_axpyl
.type blas_axpyl, @function
blas_axpyl:
## NOTE: ST(0) contains alpha

## TODO

ret

## vim: set ft=gas commentstring=##\ %s vartabstop=2,10,25,4 noexpandtab:
