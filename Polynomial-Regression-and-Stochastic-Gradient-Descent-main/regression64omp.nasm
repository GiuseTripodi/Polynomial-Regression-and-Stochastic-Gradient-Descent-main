; ---------------------------------------------------------
; Regression con istruzioni AVX a 64 bit
; ---------------------------------------------------------
; F. Angiulli
; 23/11/2017
;

;
; Software necessario per l'esecuzione:
;
;     NASM (www.nasm.us)
;     GCC (gcc.gnu.org)
;
; entrambi sono disponibili come pacchetti software 
; installabili mediante il packaging tool del sistema 
; operativo; per esempio, su Ubuntu, mediante i comandi:
;
;     sudo apt-get install nasm
;     sudo apt-get install gcc
;
; potrebbe essere necessario installare le seguenti librerie:
;
;     sudo apt-get install lib32gcc-4.8-dev (o altra versione)
;     sudo apt-get install libc6-dev-i386
;
; Per generare file oggetto:
;
;     nasm -f elf64 regression64.nasm
;

%include "sseutils64.nasm"

section .data			; Sezione contenente dati inizializzati

    align 32
    uno     dq  1.0, 1.0, 1.0, 1.0
    quattro dq  4.0

section .bss			; Sezione contenente dati non inizializzati

	alignb 32
    eta		resq	1
	condition equ 16
	condition_32 equ 32


section .text			; Sezione contenente il codice macchina

; ----------------------------------------------------------
; macro per l'allocazione dinamica della memoria
;
;	getmem	<size>,<elements>
;
; alloca un'area di memoria di <size>*<elements> bytes
; (allineata a 16 bytes) e restituisce in EAX
; l'indirizzo del primo bytes del blocco allocato
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)
;
;	fremem	<address>
;
; dealloca l'area di memoria che ha inizio dall'indirizzo
; <address> precedentemente allocata con getmem
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)

extern get_block
extern free_block

%macro	getmem	2
	mov	rdi, %1
	mov	rsi, %2
	call	get_block
%endmacro

%macro	fremem	1
	mov	rdi, %1
	call	free_block
%endmacro

global sommatoriaInterna

	;xmm0 = delta
	;rdi = Dbatchj
	;rsi = Gj
	;rdx = t
	;rcx = sommatoria

sommatoriaInterna:

    start
	vbroadcastsd  ymm8, xmm0 ; diffonde in tutto ymm0 il double di ymm0

	cmp rdx, condition
	jl aa

	sub rdx, 16
	xor rax, rax
for2:
	vmovapd ymm0, [rdi + rax * 8] ; carico Dbatch
	vmovapd ymm1, [rdi + rax * 8 + 32] ; carico Dbatch + 4
	vmovapd ymm2, [rdi + rax * 8 + 64] ; carico Dbatch + 8
	vmovapd ymm3, [rdi + rax * 8 + 96] ; carico Dbatch + 12

	vmulpd ymm0, ymm8 ; calcoliamo gjv1
	vmulpd ymm1, ymm8 ; calcoliamo gjv2
	vmulpd ymm2, ymm8 ; calcoliamo gjv3
	vmulpd ymm3, ymm8 ; calcoliamo gjv4

	vmovapd ymm4, ymm0 ; copiamo gjv1
	vmovapd ymm5, ymm1 ; copiamo gjv2
	vmovapd ymm6, ymm2 ; copiamo gjv3
	vmovapd ymm7, ymm3 ; copiamo gjv4

	vmulpd ymm0, ymm0 ;calcoliamo gjv1^2
	vmulpd ymm1, ymm1 ;calcoliamo gjv2^2
	vmulpd ymm2, ymm2 ;calcoliamo gjv3^2
	vmulpd ymm3, ymm3 ;calcoliamo gjv4^2

	vaddpd ymm0, [rsi + rax * 8] ;G[j * t + v] += powf(gjv, 2);
	vaddpd ymm1, [rsi + rax * 8 + 32] ;G[j * t + v] += powf(gjv, 2);
	vaddpd ymm2, [rsi + rax * 8 + 64 ] ;G[j * t + v] += powf(gjv, 2);
	vaddpd ymm3, [rsi + rax * 8 + 96] ;G[j * t + v] += powf(gjv, 2);

	vmovapd [rsi + rax * 8], ymm0 ; carico in memoria G
	vmovapd [rsi + rax * 8 + 32], ymm1 ; carico in memoria G
	vmovapd [rsi + rax * 8 + 64], ymm2 ; carico in memoria G
	vmovapd [rsi + rax * 8 + 96], ymm3 ; carico in memoria G

	vsqrtpd ymm0, ymm0 ;faccio la radice quadrata
	vsqrtpd ymm1, ymm1
	vsqrtpd ymm2, ymm2
	vsqrtpd ymm3, ymm3

	vdivpd ymm4 , ymm0
	vdivpd ymm5 , ymm1
	vdivpd ymm6 , ymm2
	vdivpd ymm7 , ymm3

	vaddpd  ymm4, [rcx + rax * 8 ]
	vaddpd  ymm5, [rcx + rax * 8 + 32]
	vaddpd  ymm6, [rcx + rax * 8 + 64]
	vaddpd  ymm7, [rcx + rax * 8 + 96]

	vmovapd [rcx + rax * 8 ], ymm4
	vmovapd [rcx + rax * 8  + 32], ymm5
	vmovapd [rcx + rax * 8  + 64], ymm6
	vmovapd [rcx + rax * 8  + 96], ymm7

	add rax, 16
	cmp rax, rdx
	jl for2

	add rdx, 16

aa:
	cmp rax, rdx
	jnl end

    vmovapd ymm0, [rdi + rax * 8]
	vmulpd ymm0, ymm8
	vmovapd ymm4, ymm0
	vmulpd ymm0, ymm0
	vaddpd ymm0, [rsi + rax * 8]
	vmovapd [rsi + rax * 8], ymm0
	vsqrtpd ymm0, ymm0
	vdivpd ymm4 , ymm0
	vaddpd  ymm4, [rcx + rax * 8 ]
	vmovapd [rcx + rax * 8 ], ymm4

	add rax, 4
	jmp aa

end:
	stop


;end sommatoria interna

global calcolaVettore

	;xmm0= faxt_mul (nei 64 bit meno significativi)
	;rdi= vector_src
	;rsi = vector_dest
	;rdx = lunghezza del vettore

calcolaVettore:

	push		rbp				; salva il Base Pointer
	mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
	pushaq						; salva i registri generali

	cmp rdx, condition_32
	jl ab

	sub rdx, 32
	xor rax, rax ;indice = 0
	vbroadcastsd ymm0, xmm0 ;diffonde in tutto ymm0 il double di ymm0

forcv:

	vmovapd ymm1, [rdi + rax * 8 ]
	vmovapd ymm2, [rdi + rax * 8 + 32]
	vmovapd ymm3, [rdi + rax * 8 + 64]
	vmovapd ymm4, [rdi + rax * 8 + 96]

	vmovapd ymm5, [rdi + rax * 8 + 128]
	vmovapd ymm6, [rdi + rax * 8 + 160]
	vmovapd ymm7, [rdi + rax * 8 + 192]
	vmovapd ymm8, [rdi + rax * 8 + 224]

	vmulpd ymm1, ymm1, ymm0
	vmulpd ymm2, ymm2, ymm0
	vmulpd ymm3, ymm3, ymm0
	vmulpd ymm4, ymm4, ymm0

	vmulpd ymm5, ymm5, ymm0
	vmulpd ymm6, ymm6, ymm0
	vmulpd ymm7, ymm7, ymm0
	vmulpd ymm8, ymm8, ymm0

	vaddpd ymm1, [rsi + rax * 8 ]
	vaddpd ymm2, [rsi + rax * 8 + 32 ]
	vaddpd ymm3, [rsi + rax * 8 + 64]
	vaddpd ymm4, [rsi + rax * 8 + 96]

	vaddpd ymm5, [rsi + rax * 8 + 128]
	vaddpd ymm6, [rsi + rax * 8 + 160]
	vaddpd ymm7, [rsi + rax * 8 + 192]
	vaddpd ymm8, [rsi + rax * 8 + 224]

	vmovapd [rsi + rax  * 8], ymm1 ; risultato corrente in ymm3
	vmovapd [rsi + rax  * 8 + 32], ymm2
	vmovapd [rsi + rax  * 8 + 64], ymm3
	vmovapd [rsi + rax  * 8 + 96], ymm4

	vmovapd [rsi + rax  * 8 + 128], ymm5
	vmovapd [rsi + rax  * 8 + 160], ymm6
	vmovapd [rsi + rax  * 8 + 192], ymm7
	vmovapd [rsi + rax  * 8 + 224], ymm8

	add rax, 32
	cmp rax, rdx
	jl forcv

	add rdx, 32
ab:
	cmp rax, rdx
	jnl endb

    vmovapd ymm1, [rdi + rax * 8 ]
	vmulpd ymm1, ymm1, ymm0
	vaddpd ymm1, [rsi + rax * 8 ]
	vmovapd [rsi + rax  * 8], ymm1
	add rax, 4

	jmp ab

endb:
	;stop
	popaq						; ripristina i registri generali
	mov		rsp, rbp			; ripristina lo Stack Pointer
	pop		rbp					; ripristina il Base Pointer
	ret							; torna alla funzione C chiamante

;inizio funzione prodotto scalare

global prodottoScalare

	;rdi = vector_1
	;rsi = vector_2
	;rdx = lunghezza
	;rcx = area di memeoria in cui riportare il valore del prodotto scalare

prodottoScalare:
    start

;uso un fattore di unrolling di 8, scorro quindi 32 elementi per volta

	cmp rdx, condition_32
	jl ac

	sub rdx, 32
    xor rax,rax                          ;inizializzo registro per indice corrente
    vxorpd ymm4, ymm4                       ;inizializzo vettore dove sommare i 4 (fattore unroll) registri contenentigli elementi 4 per volta

 forPR:

	vmovapd ymm0, [rdi + rax * 8]           ;metto nei registri ymm0-ymm3 4 elementi del vettore per 4 (fattore unroll) volte quindi abbiamo i primi 16 elementi
    vmovapd ymm1, [rdi + rax * 8 + 32]
   	vmovapd ymm2, [rdi + rax * 8 + 64]
    vmovapd ymm3, [rdi + rax * 8 + 96]

    vmovapd ymm5, [rdi + rax * 8 + 128]
    vmovapd ymm6, [rdi + rax * 8 + 160]
    vmovapd ymm7, [rdi + rax * 8 + 192]
    vmovapd ymm8, [rdi + rax * 8 + 224]

    vmulpd ymm0 ,ymm0,[rsi + rax * 8]             ;moltiplico i 16 elementi del primo vettore per i 16 elementi del secondo vettore
    vmulpd ymm1, ymm1, [rsi + rax * 8 + 32]
    vmulpd ymm2, ymm2, [rsi + rax * 8  + 64]
    vmulpd ymm3, ymm3, [rsi + rax * 8 + 96]

    vmulpd ymm5, ymm5, [rsi + rax * 8 + 128]
    vmulpd ymm6, ymm6, [rsi + rax * 8 + 160]
    vmulpd ymm7, ymm7, [rsi + rax * 8 + 192]
    vmulpd ymm8, ymm8, [rsi + rax * 8 + 224]

	vaddpd ymm4, ymm0;li sommo in ymm4
   	vaddpd ymm4, ymm1
    vaddpd ymm4, ymm2
    vaddpd ymm4, ymm3

    vaddpd ymm4, ymm5
    vaddpd ymm4, ymm6
    vaddpd ymm4, ymm7
    vaddpd ymm4, ymm8

	add rax, 32
    cmp rax, rdx
    jl forPR

	add rdx, 32

ac:
	cmp rax, rdx
	jnl endc

    vmovapd ymm0, [rdi + rax * 8]

	vmulpd ymm0 ,ymm0, [rsi + rax * 8]

	vaddpd ymm4, ymm0

	add rax, 4
	jmp ac

endc:

	vhaddpd ymm4, ymm4           ;halfadd x2 per sommare gli elementi del registro tra loro
	vperm2f128 ymm5, ymm4, ymm4, 00000001b
	vaddpd ymm5, ymm4

	vmovsd [rcx], xmm5

    stop

;fine funzione prodotto scalare

global inizializza
    p   equ 4
    dim equ 8
; c passata in XMM0
; lunghezza della matrice n in RSI
; x in RDI
inizializza:
        start
        SUB                     RSI, p
        VBROADCASTSD            YMM0,XMM0               ;diffondo in tutto ymm0 il double della parte bassa di xmm0
        XOR                     RCX, RCX                ;azzero l'indice
cq:     VMOVAPD                 [RDI+dim*RCX], YMM0       ;carico 4 valori di dim 8 nella matrice
        ADD                     RCX, p                  ;in ymm0 ci vanno 4 double
        CMP                     RCX, RSI
        JL                      cq
        ADD                     RSI, p
cz:     CMP                     RCX, RSI
        JNL                     fin
        VMOVSD                  [RDI+dim*RCX], XMM0       ;carico valore nella matrice
        INC                     RCX
        JMP                     cz
fin:    stop


global inserisciComb

inserisciComb:

        start
        PUSH        RBX
        PUSH        R12
        PUSH        R13
        PUSH        R14
        PUSH        R15
        VMOVAPD     YMM5, [uno]
        XOR         RAX, RAX

        XOR         R14, R14                        ;c=0
        SUB         RSI, 4                          ;RSI=numComb[h]-4    mi devo fermare quattro elementi prima
iC_c:   CMP         R14, RSI                        ;c>numComb[h]-4
        JG          iC_Rc                           ;salta al ciclo resto di ogni grado
        VMOVAPD     YMM0, YMM5                      ;accumulatore a 1.0

;		for(i=deg-1; i>=deg-h; i--) //se la combinazione e' (1,2), moltiplico ret per x[1] e poi per x[2]
        MOV         R13, RDI                        ;R13=RDI=deg
        DEC         R13                             ;R13=i=deg-1

iC_i:   CMP         R13, RDX
        JL          iC_Fi                           ;i<deg-h

;        loop vectorization
        MOV         EAX, [R8+R13*4]                ;EAX=mComb[c*deg+i] prendo un indice della combinazione
        VMOVSD      XMM2, [RCX+RAX*8]               ;prendo x[indice della combinazione] e.g. x[2]
        LEA         R12, [R8+RDI*4]                ;i+deg => raggiungo l'elemento della riga successiva

        MOV         EAX, [R12+R13*4]                ;EAX=j[indice_comb*deg+i] prendo un indice della combinazione
        VPINSRQ     XMM2, XMM2, [RCX+RAX*8], 00000001b    ;prendo x[indice della combinazione] e.g. x[2]
        LEA         R12, [R12+RDI*4]                ;i+deg => raggiungo l'elemento della riga successiva

        MOV         EAX, [R12+R13*4]                ;EAX=j[indice_comb*deg+i] prendo un indice della combinazione
        VMOVSD      XMM1, [RCX+RAX*8]               ;prendo x[indice della combinazione] e.g. x[2]
        LEA         R12, [R12+RDI*4]                ;i+deg => raggiungo l'elemento della riga successiva

        MOV         EAX, [R12+R13*4]                ;EAX=j[indice_comb*deg+i] prendo un indice della combinazione
        VPINSRQ     XMM1, XMM1, [RCX+RAX*8], 00000001b    ;prendo x[indice della combinazione] e.g. x[2]

        VPERM2F128  YMM2, YMM1, YMM2, 00000010b     ;YMM2[ XMM2 XMM1 ]
        VMULPD      YMM0, YMM2                      ;acc*=YMM2

;        fine loop vectorization
        DEC         R13                             ;i--
        JMP         iC_i


iC_fin: POP        R15
        POP        R14
        POP        R13
        POP        R12
        POP        RBX
        stop

iC_Fi:
        LEA         R8, [R12+RDI*4]                ;i+deg => raggiungo l'elemento della riga successiva
        VMOVUPD		[R9+R14*8], YMM0                ;ret[indice...indice+4]=YMM0
        ADD			R14,4			            	;c+=4
        JMP         iC_c

iC_Rc:  ADD	        RSI, 4                          ;ripristino numComb[h] che era stato diminuito

iC_RIc: CMP         R14,RSI                         ;c>=numComb[h]
        JGE         iC_fin
;		for(i=deg-1; i>=deg-h; i--)
	    MOV         R13, RDI                        ;R13=RDI=deg
        DEC         R13                             ;R13=i=deg-1
	    VMOVSD      XMM0, XMM5                      ;accumulatore a 1.0, non c'è con YMM, XMM5=YMM5?

iC_RIi: CMP         R13, RDX                        ;i<deg-h
        JL          iC_FRIi                         ;JC e' equivalente a JL
	    MOV         EAX, [R8+R13*4]                 ;EAX=j[indice_comb*deg+i] prendo un indice della combinazione
        VMULSD      XMM0, [RCX+RAX*8]               ;ac*=x[indice della combinazione]
        DEC         R13                             ;i--
        JMP         iC_RIi

iC_FRIi:
        LEA         R8, [R8+RDI*4]                  ;i+deg => raggiungo l'elemento della riga successiva
        VMOVSD      [R9+R14*8], XMM0                ;ret[indice]=XMMO[0]
        INC 		R14
        JMP         iC_RIc

global init_padd

init_padd:
        start
        VMOVUPD     YMM0, [RDI]
        SHL         RDX, 3
        DEC         RCX
        DEC         RCX

iP_c:   VMOVUPD     [RSI], YMM0
        ADD         RSI, RDX
        LOOP        iP_c

        stop