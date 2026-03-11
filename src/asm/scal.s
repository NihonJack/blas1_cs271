## scal(n, a, x) -- scale vector
## returns: x = a*x

## NOTE: Pay close attention to arguments
## 	n (INT) -> rdi
##	a (SSE) -> xmm0 / ST(0)
##	x (INT) -> rsi

## float (SSE single)
.global blas_scalf
blas_scalf:
	xor	%rax, %rax	## index = 0
	shufps	$0, %xmm0, %xmm0	## broadcast alpha from scalar x1 to packed scalar x4
		
			## FIRST perform vectorized scaling using SIMD
	movq	%rdi, %rcx	## copy n (number of floats) into rcx
	shrq	$2, %rcx	## divide rcx by 4
1:
	je	1f	## if rcx == 0, we're done
	movaps	(%rsi, %rax), %xmm1	## load 4 floats from memory into xmm1
	mulps	%xmm0, %xmm1	## multiply by alpha
	movaps	%xmm1, (%rsi, %rax)	## store the result back to memory
	leaq	16(%rax), %rax	## increment rax by 16 bytes (4x float)
	dec	%rcx
	jmp	1b
1:
			
			## NEXT perform serial scalaing using SISD
	movq	%rdi, %rcx	## copy n (number of floats) into rcx
	andq	$3, %rcx	## compute rcx mod 4
1:
	je	1f	## if rcx == 0, we're done
	movss	(%rsi, %rax), %xmm1	## load 1 float from memory into xmm1
	mulss	%xmm0, %xmm1	## multiply by alpha
	movss	%xmm1, (%rsi, %rax)	## store the result back to memory
	leaq	4(%rax), %rax	## increment rax by 8 bytes (1x float)
	dec	%rcx
	jmp	1b
1:
	movq %rsi, %rax	## return x
	ret


## double (SSE double)
.global blas_scal
.type blas_scal, @function
blas_scal:
	xor	%rax, %rax	## index = 0
	shufpd	$0, %xmm0, %xmm0	## broadcast alpha from scalar x1 to packed scalar x2

	## TODO

	movq %rsi, %rax	## return x
	ret


## long double (x87 extended floating point)
.global blas_scall
.type blas_scall, @function
blas_scall:
	## ST(0) contains alpha

	## TODO

	movq %rsi, %rax	## return x
	ret

## vim: set ft=gas commentstring=##\ %s vartabstop=2,10,25,4 noexpandtab:
