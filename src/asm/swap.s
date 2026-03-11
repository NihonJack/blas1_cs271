## swap(n, x, y) -- swap vectors

_swap:
	movq	%rdi, %rcx	## copy n (number of bytes) into rcx
	shrq	$4, %rcx	## divide rcx by 16
1:
	je	1f	## if rcx == 0, we're done
	movaps	(%rsi), %xmm0	## load 16 bytes from x into xmm0
	movaps	(%rdx), %xmm1	## load 16 bytes from y into xmm1
	
	movaps	%xmm1, (%rsi)	## store 16 bytes from xmm1 into x
	movaps	%xmm0, (%rdx)	## store 16 bytes from xmm0 into y
	
	leaq	16(%rsi), %rsi	## increment rdx by 16 bytes
	leaq	16(%rdx), %rdx	## increment rdx by 16 bytes
	dec	%rcx
	jmp	1b
1:

	test	$0b1000, %rdi
	je 1f
	movq	(%rsi), %rax
	movq	(%rdx), %rcx

	movq	%rax, (%rdx)
	movq	%rcx, (%rsi)
	leaq	8(%rsi), %rsi
	leaq	8(%rdx), %rdx
1:

	test $0b100, %rdi
	je 1f
	movl	(%rsi), %eax
	movl	(%rdx), %ecx

	movl	%eax, (%rdx)
	movl	%ecx, (%rsi)
	leaq	4(%rsi), %rsi
	leaq	4(%rdx), %rdx
1:

	test $0b10, %rdi
	je 1f
	movw	(%rsi), %ax
	movw	(%rdx), %cx

	movw	%ax, (%rdx)
	movw	%cx, (%rsi)
	leaq	2(%rsi), %rsi
	leaq	2(%rdx), %rdx
1:

	test $0b1, %rdi
	je 1f
	movb	(%rsi), %al
	movb	(%rdx), %cl

	movb	%al, (%rdx)
	movb	%cl, (%rsi)
1:
	ret

## float (SSE single)
.global blas_swapf
.type blas_swapf, @function
blas_swapf:
	shlq $2, %rdi
	jmp _swap

## double (SSE double)
.global blas_swap
.type blas_swap, @function
blas_swap:
	shlq $3, %rdi
	jmp _swap


## long double (x87 extended floating point)
.global blas_swapl
.type blas_swapl, @function
blas_swapl:
	shlq $4, %rdi
	jmp _swap

## vim: set ft=gas commentstring=##\ %s vartabstop=2,10,25,4 noexpandtab:
