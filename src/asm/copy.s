## copy(n, x, y) -- copy vector
## y = x

## float (SSE single)
.global blas_copyf
.type blas_copyf, @function
blas_copyf:
	shlq $2, %rdi
	movq %rdi, %rax
	movq %rdx, %rdi
	movq %rax, %rdx
	jmp memcpy@plt

## double (SSE double)
.global blas_copy
.type blas_copy, @function
blas_copy:
	shlq $3, %rdi
	movq %rdi, %rax
	movq %rdx, %rdi
	movq %rax, %rdx
	jmp memcpy@plt


## long double (x87 extended floating point)
.global blas_copyl
.type blas_copyl, @function
blas_copyl:
	shlq $4, %rdi
	movq %rdi, %rax
	movq %rdx, %rdi
	movq %rax, %rdx
	jmp memcpy@plt

## vim: set ft=gas commentstring=##\ %s vartabstop=2,10,25,4 noexpandtab:
