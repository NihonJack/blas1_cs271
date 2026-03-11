## asum(n, x) -- 1-norm
## returns: |x|

## NOTE: This implementation has been completed for you as an example

## float (SSE single)
.global blas_asumf
.type blas_asumf, @function
blas_asumf:
	xor	%rax, %rax	## index = 0 
	xorps	%xmm0, %xmm0	## clear accumulator

			## xmm2 will be used to mask sign bits
	pcmpeqd	%xmm2, %xmm2	## set xmm2 all 1s (0xFFFFFFFF x4)
	psrld	$1, %xmm2	## shift left 1    (0x7FFFFFFF x4)

			## Vectorized summation using SIMD
	movq	%rdi, %rcx	## copy n (number of floats) into rcx
	shrq	$2, %rcx	## divide rcx by 4
1:
	je	1f	## if rcx == 0, we're done
	movaps	(%rsi, %rax), %xmm1	## load 4 floats from memory into xmm1
	andps	%xmm2, %xmm1	## mask the sign bits (absolute value)
	addps	%xmm1, %xmm0	## accumulate into xmm0
	leaq	16(%rax), %rax	## increment rax by 16 bytes (4x float)
	dec	%rcx
	jmp	1b
1:
			## Horizontal reduction of xmm0
	movhlps	%xmm0, %xmm1	## xmm1[0:1] = xmm0[2:3]
	addps	%xmm1, %xmm0	## xmm0[0:1] += xmm1[0:1]
	pshufd	$1, %xmm0, %xmm1	## xmm1[0] = xmm0[1]
	addss	%xmm1, %xmm0	## xmm0[0] += xmm1[0]
			## xmm0 now contains a single scalar value

			## Serial summation using SISD
	movq	%rdi, %rcx	## copy n (number of floats) into rcx
	andq	$3, %rcx	## compute rcx mod 4
1:
	je	1f	## if rcx == 0, we're done
	movss	(%rsi, %rax), %xmm1	## load 1 float from memory into xmm1
	andps	%xmm2, %xmm1	## mask the sign bit (absolute value)
	addss	%xmm1, %xmm0	## accumulate into xmm0
	leaq	4(%rax), %rax	## increment rax by 8 bytes (1x float)
	dec	%rcx
	jmp	1b
1:
ret


## double (SSE double)
.global blas_asum
.type blas_asum, @function
blas_asum:
	xor	%rax, %rax	## index = 0 
	xorps	%xmm0, %xmm0	## clear accumulator

			## xmm2 will be used to mask sign bits
	pcmpeqd	%xmm2, %xmm2	## set xmm2 all 1s (0xFF...FF x2)
	psrlq	$1, %xmm2	## shift left 1    (0x7F...FF x2)
			
			## Vectorized summation using SIMD
	movq	%rdi, %rcx	## copy n (number of doubles) into rcx
	shrq	$1, %rcx	## divide rcx by 2
1:
	je	1f	## if rcx == 0, we're done
	movapd	(%rsi, %rax), %xmm1	## load 2 doubles from memory into xmm1
	andpd	%xmm2, %xmm1	## mask the sign bits (absolute value)
	addpd	%xmm1, %xmm0	## accumulate into xmm0
	leaq	16(%rax), %rax	## increment rax by 16 bytes (2x double)
	dec	%rcx
	jmp	1b
1:
			## Horizontal reduction of xmm0
	movhlps	%xmm0, %xmm1 	## xmm1[0] = xmm0[1]
	addsd	%xmm1, %xmm0 	## xmm0[0] += xmm1[0]
			## xmm0 now contains a single scalar value
			
			## Serial summation using SISD
	movq	%rdi, %rcx	## copy n (number of floats) into rcx
	andq	$1, %rcx	## compute rcx mod 2
	je	1f	## if rcx == 0, we're done
	movsd	(%rsi, %rax), %xmm1	## load 1 double from memory into xmm1
	andpd	%xmm2, %xmm1	## mask the sign bit (absolute value)
	addsd	%xmm1, %xmm0	## accumulate into xmm0
1:
ret


## long double (x87 extended floating point)
.global blas_asuml
.type blas_asuml, @function
blas_asuml:
	xor	%rax, %rax	## index = 0

	fldz		## ST0 = 0.0 (accumulator)
	test	%rdi, %rdi
1:
	je	1f
	fldt	(%rsi, %rax)	## load x[i]
	fabs
	faddp	%st(1)	## Add ST1 to ST0
	leaq	16(%rax), %rax	## index += 16
	dec	%rdi
	jmp	1b
1:
ret

## vim: set ft=gas commentstring=##\ %s vartabstop=2,10,25,4 noexpandtab:
