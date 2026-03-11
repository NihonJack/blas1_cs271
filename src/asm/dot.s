## dot(n, x, y) -- dot product
## returns: x \cdot y


## float (SSE single)
.global blas_dotf
blas_dotf:

## TODO

ret


## double (SSE double)
.global blas_dot
.type blas_dot, @function
blas_dot:

## TODO

ret


## long double (x87 extended floating point)
.global blas_dotl
.type blas_dotl, @function
blas_dotl:

## TODO

ret

## vim: set ft=gas commentstring=##\ %s vartabstop=2,10,25,4 noexpandtab:
