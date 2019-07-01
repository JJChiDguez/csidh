
.intel_syntax noprefix

.section .rodata

.set pbits, 511
.global p
p:
    .quad 0x1b81b90533c6c87b, 0xc2721bf457aca835, 0x516730cc1f0b4f25, 0xa7aac6c567f35507
    .quad 0x5afbfcc69322c9cd, 0xb42d083aedc88c42, 0xfc8ab0d15e3e4c4a, 0x65b48e8f740f89bf


.global R_mod_p
R_mod_p: /* 2^512 mod p */
    .quad 0xc8fc8df598726f0a, 0x7b1bc81750a6af95, 0x5d319e67c1e961b4, 0xb0aa7275301955f1
    .quad 0x4a080672d9ba6c64, 0x97a5ef8a246ee77b, 0x06ea9e5d4383676a, 0x3496e2e117e0ec80


/* (2^512)^2 mod p */
.global R_squared_mod_p
R_squared_mod_p:
    .quad 0x36905b572ffc1724, 0x67086f4525f1f27d, 0x4faf3fbfd22370ca, 0x192ea214bcc584b1
    .quad 0x5dae03ee2f5de3d0, 0x1e9248731776b371, 0xad5f166e20e4f52d, 0x4ed759aea6f3917e

/* -p^-1 mod 2^64 */
.inv_min_p_mod_r:
    .quad 0x66c1301f632e294d


.section .text

/* Integer addition */
.u512_add:
    mov rax, [rsi +  0]
    add rax, [rdx +  0]
    mov [rdi +  0], rax
    .set k, 1
    .rept 7
        mov rax, [rsi + 8*k]
        adc rax, [rdx + 8*k]
        mov [rdi + 8*k], rax
        .set k, k+1
    .endr
    setc al
    movzx rax, al
    ret

/* Integer substraction */
.u512_sub:
    mov rax, [rsi +  0]
    sub rax, [rdx +  0]
    mov [rdi +  0], rax
    .set k, 1
    .rept 7
        mov rax, [rsi + 8*k]
        sbb rax, [rdx + 8*k]
        mov [rdi + 8*k], rax
        .set k, k+1
    .endr
    setc al
    movzx rax, al
    ret

.global fp_copy
fp_copy:
    cld
    mov rcx, 8
    rep movsq
    ret

.global fp_cswap
fp_cswap:
    movzx rax, dl
    neg rax
    .set k, 0
    .rept 8
        mov rcx, [rdi + 8*k]
        mov rdx, [rsi + 8*k]

        mov r8, rcx
        xor r8, rdx
        and r8, rax

        xor rcx, r8
        xor rdx, r8

        mov [rdi + 8*k], rcx
        mov [rsi + 8*k], rdx

        .set k, k+1
    .endr
    ret

.reduce_once:
    push rbp
    mov rbp, rdi

    mov rdi, [rbp +  0]
    sub rdi, [rip + p +  0]
    mov rsi, [rbp +  8]
    sbb rsi, [rip + p +  8]
    mov rdx, [rbp + 16]
    sbb rdx, [rip + p + 16]
    mov rcx, [rbp + 24]
    sbb rcx, [rip + p + 24]
    mov r8,  [rbp + 32]
    sbb r8,  [rip + p + 32]
    mov r9,  [rbp + 40]
    sbb r9,  [rip + p + 40]
    mov r10, [rbp + 48]
    sbb r10, [rip + p + 48]
    mov r11, [rbp + 56]
    sbb r11, [rip + p + 56]

    setnc al
    movzx rax, al
    neg rax

.macro cswap2, r, m
    xor \r, \m
    and \r, rax
    xor \m, \r
.endm

    cswap2 rdi, [rbp +  0]
    cswap2 rsi, [rbp +  8]
    cswap2 rdx, [rbp + 16]
    cswap2 rcx, [rbp + 24]
    cswap2 r8,  [rbp + 32]
    cswap2 r9,  [rbp + 40]
    cswap2 r10, [rbp + 48]
    cswap2 r11, [rbp + 56]

    pop rbp
    ret

.global fp_add
fp_add:
    push rdi
    call .u512_add
    pop rdi
    jmp .reduce_once

.global fp_sub
fp_sub:
    push rdi
    call .u512_sub
    pop rdi
    xor rsi, rsi
    xor rdx, rdx
    xor rcx, rcx
    xor r8, r8
    xor r9, r9
    xor r10, r10
    xor r11, r11
    test rax, rax
    cmovnz rax, [rip + p +  0]
    cmovnz rsi, [rip + p +  8]
    cmovnz rdx, [rip + p + 16]
    cmovnz rcx, [rip + p + 24]
    cmovnz r8,  [rip + p + 32]
    cmovnz r9,  [rip + p + 40]
    cmovnz r10, [rip + p + 48]
    cmovnz r11, [rip + p + 56]
    add [rdi +  0], rax
    adc [rdi +  8], rsi
    adc [rdi + 16], rdx
    adc [rdi + 24], rcx
    adc [rdi + 32],  r8
    adc [rdi + 40],  r9
    adc [rdi + 48], r10
    adc [rdi + 56], r11
    ret

/* Montgomery arithmetic */

.global fp_mul
fp_mul:
    push rbp
    push rbx
    push r12
    push r13
    push r14
    push r15

    push rdi

    mov rdi, rsi
    mov rsi, rdx

    xor r8,  r8
    xor r9,  r9
    xor r10, r10
    xor r11, r11
    xor r12, r12
    xor r13, r13
    xor r14, r14
    xor r15, r15
    xor rbp, rbp

    /* flags are already cleared */

.macro MULSTEP, k, r0, r1, r2, r3, r4, r5, r6, r7, r8

    mov rdx, [rsi +  0]
    mulx rcx, rdx, [rdi + 8*\k]
    add rdx, \r0
    mulx rcx, rdx, [rip + .inv_min_p_mod_r]

    xor rax, rax /* clear flags */

    mulx rbx, rax, [rip + p +  0]
    adox \r0, rax

    mulx rcx, rax, [rip + p +  8]
    adcx \r1, rbx
    adox \r1, rax

    mulx rbx, rax, [rip + p + 16]
    adcx \r2, rcx
    adox \r2, rax

    mulx rcx, rax, [rip + p + 24]
    adcx \r3, rbx
    adox \r3, rax

    mulx rbx, rax, [rip + p + 32]
    adcx \r4, rcx
    adox \r4, rax

    mulx rcx, rax, [rip + p + 40]
    adcx \r5, rbx
    adox \r5, rax

    mulx rbx, rax, [rip + p + 48]
    adcx \r6, rcx
    adox \r6, rax

    mulx rcx, rax, [rip + p + 56]
    adcx \r7, rbx
    adox \r7, rax

    mov rax, 0
    adcx \r8, rcx
    adox \r8, rax


    mov rdx, [rdi + 8*\k]

    xor rax, rax /* clear flags */

    mulx rbx, rax, [rsi +  0]
    adox \r0, rax

    mulx rcx, rax, [rsi +  8]
    adcx \r1, rbx
    adox \r1, rax

    mulx rbx, rax, [rsi + 16]
    adcx \r2, rcx
    adox \r2, rax

    mulx rcx, rax, [rsi + 24]
    adcx \r3, rbx
    adox \r3, rax

    mulx rbx, rax, [rsi + 32]
    adcx \r4, rcx
    adox \r4, rax

    mulx rcx, rax, [rsi + 40]
    adcx \r5, rbx
    adox \r5, rax

    mulx rbx, rax, [rsi + 48]
    adcx \r6, rcx
    adox \r6, rax

    mulx rcx, rax, [rsi + 56]
    adcx \r7, rbx
    adox \r7, rax

    mov rax, 0
    adcx \r8, rcx
    adox \r8, rax

.endm

    MULSTEP 0, r8,  r9,  r10, r11, r12, r13, r14, r15, rbp
    MULSTEP 1, r9,  r10, r11, r12, r13, r14, r15, rbp, r8
    MULSTEP 2, r10, r11, r12, r13, r14, r15, rbp, r8,  r9
    MULSTEP 3, r11, r12, r13, r14, r15, rbp, r8,  r9,  r10
    MULSTEP 4, r12, r13, r14, r15, rbp, r8,  r9,  r10, r11
    MULSTEP 5, r13, r14, r15, rbp, r8,  r9,  r10, r11, r12
    MULSTEP 6, r14, r15, rbp, r8,  r9,  r10, r11, r12, r13
    MULSTEP 7, r15, rbp, r8,  r9,  r10, r11, r12, r13, r14

    pop rdi

    mov [rdi +  0], rbp
    mov [rdi +  8], r8
    mov [rdi + 16], r9
    mov [rdi + 24], r10
    mov [rdi + 32], r11
    mov [rdi + 40], r12
    mov [rdi + 48], r13
    mov [rdi + 56], r14

    pop r15
    pop r14
    pop r13
    pop r12
    pop rbx
    pop rbp
    jmp .reduce_once

.fp_mul2:
    mov rdx, rdi
    jmp fp_mul

.global fp_sqr
fp_sqr:
    /* TODO implement optimized Montgomery squaring */
    mov rdx, rsi
    jmp fp_mul

.fp_sq1:
    mov rsi, rdi
    jmp fp_sqr

/* (obviously) not constant time in the exponent! */
.fp_pow:
    push rbx
    mov rbx, rsi
    push r12
    push r13
    push rdi
    sub rsp, 64

    mov rsi, rdi
    mov rdi, rsp
    call fp_copy

    mov rdi, [rsp + 64]
    lea rsi, [rip + R_mod_p]
    call fp_copy

.macro POWSTEP, k
        mov r13, [rbx + 8*\k]
        xor r12, r12

        0:
        test r13, 1
        jz 1f

        mov rdi, [rsp + 64]
        mov rsi, rsp
        call .fp_mul2

        1:
        mov rdi, rsp
        call .fp_sq1

        shr r13

        inc r12
        test r12, 64
        jz 0b
.endm

    POWSTEP 0
    POWSTEP 1
    POWSTEP 2
    POWSTEP 3
    POWSTEP 4
    POWSTEP 5
    POWSTEP 6
    POWSTEP 7

    add rsp, 64+8
    pop r13
    pop r12
    pop rbx
    ret

.section .rodata
.p_minus_2:
    .quad 0x1b81b90533c6c879, 0xc2721bf457aca835, 0x516730cc1f0b4f25, 0xa7aac6c567f35507
    .quad 0x5afbfcc69322c9cd, 0xb42d083aedc88c42, 0xfc8ab0d15e3e4c4a, 0x65b48e8f740f89bf

.section .text

/* TODO use a better addition chain? */
.global fp_inv
fp_inv:
    lea rsi, [rip + .p_minus_2]
    jmp .fp_pow

.section .rodata
.global p_minus_1_halves
p_minus_1_halves:
    .quad 0x8dc0dc8299e3643d, 0xe1390dfa2bd6541a, 0xa8b398660f85a792, 0xd3d56362b3f9aa83
    .quad 0x2d7dfe63499164e6, 0x5a16841d76e44621, 0xfe455868af1f2625, 0x32da4747ba07c4df

.section .text

/* TODO use a better addition chain? */
.global fp_issquare
fp_issquare:
    push rdi
    lea rsi, [rip + p_minus_1_halves]
    call .fp_pow
    pop rdi

    xor rax, rax
    .set k, 0
    .rept 8
        mov rsi, [rdi + 8*k]
        xor rsi, [rip + R_mod_p + 8*k]
        or rax, rsi
        .set k, k+1
    .endr
    test rax, rax
    setz al
    movzx rax, al
    ret


/* not constant time (but this shouldn't leak anything of importance) */
.global fp_random
fp_random:

    push rdi
    mov rsi, 64
    call randombytes
    pop rdi
    mov rax, 1
    shl rax, (pbits % 64)
    dec rax
    and [rdi + 56], rax

    .set k, 7
    .rept 8
        mov rax, [rip + p + 8*k]
        cmp [rdi + 8*k], rax
        jge fp_random
        jl 0f
        .set k, k-1
    .endr
    0:
    ret

