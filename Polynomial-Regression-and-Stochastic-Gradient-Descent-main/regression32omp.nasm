; ---------------------------------------------------------
; Regressione con istruzioni SSE a 32 bit
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
;     nasm -f elf32 regression32.nasm 
;
%include "sseutils32.nasm"

section .data			; Sezione contenente dati inizializzati
    align 16
    uno     dd  1.0, 1.0, 1.0, 1.0
    quattro dd  4.0
	condition equ 16
	condition_24 equ 24 ; condizione di almeno 24 elementi nel vettore


section .bss			; Sezione contenente dati non inizializzati
	alignb 16
	eta		resd		1
	alignb 16
	delta resd 4



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
	mov	eax, %1
	push	eax
	mov	eax, %2
	push	eax
	call	get_block
	add	esp, 8
%endmacro

%macro	fremem	1
	push	%1
	call	free_block
	add	esp, 4
%endmacro

; ------------------------------------------------------------
; Funzioni
; ------------------------------------------------------------


global sommatoriaInterna

	delta_param equ 8
	Dbatch equ 12
	G equ 16
	t equ 20
	sommatoria equ 24



sommatoriaInterna:

	start

	mov ebx, [ebp + Dbatch] ; ebx = Dbatch
	mov ecx, [ebp + G]; ecx = puntatore G
	mov edi, [ebp + t]; edi = t
	mov edx, [ebp + sommatoria]

	cmp edi, condition
	jl rt

	movss xmm0, [ebp  + delta_param];xmm0 = delta
	shufps  xmm0, xmm0, 0 ;diffonde in tutto xmm0 il double di xmm0
	movaps [delta], xmm0

	sub edi, 16
	xor eax, eax
for2:
	movaps xmm0, [ebx + eax * 4] ; carico Dbatch
	movaps xmm1, [ebx + eax * 4 + 16] ; carico Dbatch + 4
	movaps xmm2, [ebx + eax * 4 + 32] ; carico Dbatch + 8
	movaps xmm3, [ebx + eax * 4 + 48] ; carico Dbatch + 12

	mulps xmm0, [delta] ; calcoliamo gjv1
	mulps xmm1, [delta] ; calcoliamo gjv2
	mulps xmm2, [delta] ; calcoliamo gjv3
	mulps xmm3, [delta] ; calcoliamo gjv4

	movaps xmm4, xmm0 ; copiamo gjv1
	movaps xmm5, xmm1 ; copiamo gjv2
	movaps xmm6, xmm2 ; copiamo gjv3
	movaps xmm7, xmm3 ; copiamo gjv4

	mulps xmm0, xmm0 ;calcoliamo gjv1^2
	mulps xmm1, xmm1 ;calcoliamo gjv2^2
	mulps xmm2, xmm2 ;calcoliamo gjv3^2
	mulps xmm3, xmm3 ;calcoliamo gjv4^2

	addps xmm0, [ecx + eax * 4] ;G[j * t + v] += powf(gjv, 2);
	addps xmm1, [ecx + eax * 4 + 16] ;G[j * t + v] += powf(gjv, 2);
	addps xmm2, [ecx + eax * 4 + 32 ] ;G[j * t + v] += powf(gjv, 2);
	addps xmm3, [ecx + eax * 4 + 48] ;G[j * t + v] += powf(gjv, 2);

	movaps [ecx + eax * 4], xmm0 ; carico in memoria G
	movaps [ecx + eax * 4 + 16], xmm1 ; carico in memoria G
	movaps [ecx + eax * 4 + 32], xmm2 ; carico in memoria G
	movaps [ecx + eax * 4 + 48], xmm3 ; carico in memoria G


	sqrtps xmm0, xmm0 ;faccio la radice quadrata
	sqrtps xmm1, xmm1
	sqrtps xmm2, xmm2
	sqrtps xmm3, xmm3

	divps xmm4 , xmm0
	divps xmm5 , xmm1
	divps xmm6 , xmm2
	divps xmm7 , xmm3

	addps  xmm4, [edx + eax * 4]
	addps  xmm5, [edx + eax * 4 + 16]
	addps  xmm6, [edx + eax * 4 + 32]
	addps  xmm7, [edx + eax * 4+ 48]

	movaps [edx + eax * 4], xmm4
	movaps [edx + eax * 4 + 16], xmm5
	movaps [edx + eax * 4 + 32], xmm6
	movaps [edx + eax * 4 + 48], xmm7

	add eax, 16
	cmp eax, edi
	jl for2

	add edi, 16
	movaps xmm7, [delta]
rt:

	cmp eax, edi
	jnl enda

	movaps xmm0, [ebx + eax * 4] ; carico Dbatch
	mulps xmm0, xmm7 ; calcoliamo gjv1
	movaps xmm4, xmm0 ; copiamo gjv1
	mulps xmm0, xmm0 ;calcoliamo gjv1^2
	addps xmm0, [ecx + eax * 4] ;G[j * t + v] += powf(gjv, 2);
	movaps [ecx + eax * 4], xmm0 ; carico in memoria G
	sqrtps xmm0, xmm0 ;faccio la radice quadrata
	divps xmm4 , xmm0
	addps  xmm4, [edx + eax * 4]
	movaps [edx + eax * 4], xmm4

	add eax, 4
	jmp rt


enda:

stop
;fine funzione sommatoria interna

global calcolaVettore

	fact_mul equ 8
	vector_src equ 12
	vector_dest equ 16
	len equ 20

calcolaVettore:
	start

	mov ebx, [ebp + vector_src] ; ebx = vector_src
	mov ecx, [ebp + vector_dest]; ecx = puntatore a vector_dest
	mov edi, [ebp + len]; edi = len

	cmp edi, condition_24 ; verifico che ci siano alemeno 24 elementi nel vettore
	jl rs

	movss xmm0, [ebp  + fact_mul] ;xmm0 = fact_mul; fattore moltiplicativo

	xor eax, eax ;indice = 0
	shufps  xmm0, xmm0, 0 ;diffonde in tutto ymm0 il float di xmm0

    ;il fattore moltiplicativo è 6, quindi nell'unrolling scorro 24 elementi per volta
	sub edi, 24

for1:

	movaps xmm1, [ebx + eax * 4]    ; carico il vect_src
	movaps xmm2, [ebx + eax * 4 + 16]
	movaps xmm3, [ebx + eax * 4 + 32]
	movaps xmm4, [ebx + eax * 4 + 48]
	movaps xmm5, [ebx + eax * 4 + 64]
	movaps xmm6, [ebx + eax * 4 + 80]


	mulps xmm1, xmm0 ; moltiplica il facr_mul per vecto_src[j :j + 4]
	mulps xmm2, xmm0
	mulps xmm3, xmm0
	mulps xmm4, xmm0
	mulps xmm5, xmm0
	mulps xmm6, xmm0

	addps xmm1, [ecx + eax * 4] ;valore di sommatoria
	addps xmm2, [ecx + eax * 4 + 16]
	addps xmm3, [ecx + eax * 4 + 32]
	addps xmm4, [ecx + eax * 4 + 48]
	addps xmm5, [ecx + eax * 4 + 64]
	addps xmm6, [ecx + eax * 4 + 80]

	movaps [ecx + eax  * 4], xmm1 ; risultato corrente in ymm3
	movaps [ecx + eax  * 4 + 16], xmm2
	movaps [ecx + eax  * 4 + 32], xmm3
	movaps [ecx + eax  * 4 + 48], xmm4
	movaps [ecx + eax  * 4 + 64], xmm5
	movaps [ecx + eax  * 4 + 80], xmm6

	add eax, 24;ci muoviamo di 16 elementi per volta
	cmp eax, edi
	jl for1

	add edi , 24
rs:

	cmp eax, edi
	jnl endb

	movaps xmm1, [ebx + eax * 4 ]    ; carico il vect_src
	mulps xmm1, xmm0 ; moltiplica il facr_mul per vecto_src[j :j + 4]
	addps xmm1, [ecx + eax * 4 ] ;valore di sommatoria
	movaps [ecx + eax  * 4], xmm1 ; risultato corrente in ymm3

	add eax, 4;ci muoviamo di 16 elementi per volta
	jmp rs

endb:
	stop

global prodottoScalare

	;predo i parametri della funzione

	ps_A equ 8
	ps_B equ 12
	ps_dim equ 16
	ps_value_ret equ 20

prodottoScalare:
	;sequenza di ingresso nella funzione

	push ebp
	mov ebp, esp
	push ebx
	push esi
	push edi

	;prendo i parametri

	mov ebx, [ebp + ps_A] ;ebx = puntatore al vettore A
	mov ecx, [ebp + ps_B] ;ecx = puntatore al vettore B
	mov edi, [ebp + ps_dim ]; edi = dimensione del vettore
	mov esi, [ebp + ps_value_ret]; xmm3 = nei primi 32 vit puntatore al valore di ritorno

	cmp edi, condition_24
	jl rp

	xor eax, eax				; i = 0 azzero il contatore
	xorps xmm4, xmm4

	sub edi, 24 ;leggo 24 elementi del vettore per volta
for:
	movaps xmm0, [ebx + eax * 4] ;prendo i primi  4 elementi dal primo vettore
	movaps xmm1, [ebx + eax * 4 + 16] ;prendo i secondi 4 elementi dal primo vettore
	movaps xmm2, [ebx + eax * 4 + 32] ;prendo i terzi 4 elementi dal primo vettore
	movaps xmm3, [ebx + eax * 4 + 48] ;prendo i quarti 4 elementi dal primo vettore
	movaps xmm5, [ebx + eax * 4 + 64]
	movaps xmm6, [ebx + eax * 4 + 80]

	mulps xmm0,[ecx + eax * 4 ] ;prendo 4 elementi dal secondo vettore e moltiplico
	mulps xmm1,[ecx + eax * 4 + 16] ;prendo 4 elementi dal secondo vettore e moltiplico
	mulps xmm2,[ecx + eax * 4 + 32] ;prendo 4 elementi dal secondo vettore e moltiplico
	mulps xmm3,[ecx + eax * 4 + 48] ;prendo 4 elementi dal secondo vettore e moltiplico
	mulps xmm5,[ecx + eax * 4 + 64]
	mulps xmm6,[ecx + eax * 4 + 80]

	addps xmm2, xmm0 ;metto il risultato corrente in xmm4
	addps xmm3, xmm1 ;metto il risultato corrente in xmm4
	addps xmm4, xmm2 ;metto il risultato corrente in xmm4
	addps xmm4, xmm3 ;metto il risultato corrente in xmm4
	addps xmm4, xmm5 ;metto il risultato corrente in xmm4
	addps xmm4, xmm6 ;metto il risultato corrente in xmm4

	add eax, 24 ; mi muovo di 16 elementi per volta se il fattore di unrolling e 4

	cmp eax, edi
	jl for

	add edi, 24
rp:

	cmp eax, edi
	jnl endd

	movaps xmm0, [ebx + eax * 4] ;prendo i primi  4 elementi dal primo vettore
	mulps xmm0,[ecx + eax * 4 ] ;prendo 4 elementi dal secondo vettore e moltiplico

	addps xmm4, xmm0 ;metto il risultato corrente in xmm4

	add eax, 4 ; mi muovo di 16 elementi per volta se il fattore di unrolling e 4


	jmp rp

endd:

	;compatto il risultato
	haddps xmm4, xmm4
	haddps xmm4, xmm4

	movss [esi], xmm4

	;sequenza di uscita dalla funzione

	pop edi
	pop esi
	pop ebx
	mov esp, ebp
	pop ebp
	ret

global inizializza
    x   equ 8   ;matrice
    n   equ 12  ;lunghezza della matrice
    con equ 16
    un  equ 1

inizializza:
        start
        MOV             EAX, [EBP+x]
        MOV             EDI, [EBP+n]
        SUB             EDI, 4*un
        MOVSS           XMM0, [EBP+con]
        SHUFPS          XMM0, XMM0, 00000000b
        XOR             ESI, ESI
    cq: MOVAPS          [EAX+4*ESI], XMM0       ;carico valore nella matrice
        ADD             ESI, 4*un
        CMP             ESI, EDI
        JL              cq
        ADD             EDI, 4*un
    cz: CMP             ESI, EDI
        JNL             fin
        MOVSS           [EAX+4*ESI], XMM0   ;carico valore nella matrice
        INC             ESI
        JMP             cz
    fin:stop

global inserisciComb

    iC_deg          equ 8
    iC_nCH          equ 12
    iC_degh         equ 16
    iC_mComb        equ 20
    iC_xast         equ 24
    iC_x            equ 28

    shiftLC3        equ 11001001b
    shiftLCZ3       equ 11001011b
    shiftRC3        equ 11010010b
    shift2C         equ 01001110b
    blend34         equ 00001100b
    shiftRC         equ 10010011b
    shiftLC         equ 00111001b
    shift2LC         equ 11110001b

inserisciComb:

        start
        MOVAPS      XMM5, [uno]
        CVTPI2PS    XMM7, [EBP+iC_nCH]                ;XMM7[ iC_nCH iC_degh - - ]
        XORPS       XMM6, XMM6
        MOV         EDX, [EBP+iC_deg]
        DEC         EDX
        CVTSI2SS    XMM4, EDX                       ;XMM4[0]=deg-1
        INC         EDX
        MOV         EBX, [EBP+iC_xast]
        MOV         EAX, [EBP+iC_x]
        MOV         ECX, [EBP+iC_mComb]
        MOVSS       XMM3, [quattro]

;       for(c=0; c<numComb[h]; c++) //c � J, per ogni J appartenente a [d]^h
        SUBSS       XMM7, XMM3                      ;XMM7[0]=numComb[h]-4    mi devo fermare quattro elementi prima
iC_c:   COMISS      XMM6, XMM7                      ;c>numComb[h]-4
        JA          iC_Rc                           ;salta al ciclo resto di ogni grado JA è equivalente a JG
        MOVAPS      XMM0, XMM5                      ;accumulatore a 1.0

;		for(i=deg-1; i>=deg-h; i--) //se la combinazione � (1,2), moltiplico ret per x[1] e poi per x[2]
        SHUFPS      XMM7, XMM7, shift2LC            ;XMM7[0]=deg-h
        SHUFPS      XMM6, XMM6, shift2LC            ;XMM6[0]<=>i
        MOVSS       XMM6, XMM4                      ;XMM6[0]=deg-1 <=> i=deg-1

iC_i:   COMISS      XMM6, XMM7                      ;i<deg-h
        JC          iC_Fi                           ;JC e' equivalente a JL

        CVTTSS2SI   EDI, XMM6                       ;MOVSS EDI, XMM6 <=> EDI=i
;        loop vectorization
        MOV         ESI, [ECX+EDI*4]                ;ESI=j[indice_comb*deg+i] prendo un indice della combinazione
        INSERTPS    XMM1, [EAX+ESI*4], 00000000b    ;prendo x[indice della combinazione] e.g. x[2]
        ADD         EDI, EDX                        ;i+deg => raggiungo l'elemento della riga successiva
        MOV         ESI, [ECX+EDI*4]                ;ESI=j[indice_comb*deg+i]
        INSERTPS    XMM1, [EAX+ESI*4], 00010000b
        ADD         EDI, EDX
        MOV         ESI, [ECX+EDI*4]                ;ESI=j[indice_comb*deg+i]
        INSERTPS    XMM1, [EAX+ESI*4], 00100000b
        ADD         EDI, EDX
        MOV         ESI, [ECX+EDI*4]                ;ESI=j[indice_comb*deg+i]
        INSERTPS    XMM1, [EAX+ESI*4], 00110000b
        MULPS       XMM0, XMM1                      ;accumulatore *= colonne correnti
;        fine loop vectorization
        SUBSS       XMM6, XMM5                      ;XMM6[0]-- <=> i--
        JMP         iC_i


iC_fin: stop

iC_Fi:  SHUFPS      XMM7, XMM7, shift2LC            ;XMM7[0]=numComb[h]
        SHUFPS      XMM6, XMM6, shift2LC            ;XMM6[0]=c
        ADDSS       XMM6, XMM3                      ;c+=4
        MOV         EDI, EDX                        ;EDI=deg
        SHL         EDI, 4                          ;EDI*16
        ADD         ECX, EDI                        ;ECX+deg*16, mi sposto di 4 righe su mComb
        MOVUPS      [EBX], XMM0               ;ret[indice...indice+4]=XMM0
        ADD         EBX, 16
        JMP         iC_c


iC_Rc:  ADDSS       XMM7, XMM3                      ;ripristino numComb[h] che era stato diminuito
iC_RIc: COMISS      XMM6, XMM7                      ;c>=numComb[h]
        JNC         iC_fin

;		for(i=deg-1; i>=deg-h; i--)
        SHUFPS      XMM7, XMM7, shift2LC           ;XMM7[0]=deg-h
        SHUFPS      XMM6, XMM6, shift2LC           ;XMM6[0]<=>i
        MOVSS       XMM6, XMM4                      ;XMM6[0]=deg-1 <=> i=deg-1
        MOVSS       XMM0, XMM5                      ;accumulatore a 1.0
iC_RIi: COMISS      XMM6, XMM7                      ;i<deg-h
        JC          iC_FRIi                         ;JC e' equivalente a JL
        CVTTSS2SI   EDI, XMM6                       ;MOVSS EDI, XMM6 <=> EDI=i
        MOV         ESI, [ECX+EDI*4]                ;ESI=j[indice_comb*deg+i] prendo un indice della combinazione
        MULSS       XMM0, [EAX+ESI*4]               ;prendo x[indice della combinazione]
        SUBSS       XMM6, XMM5                      ;XMM6[0]-- <=> i--
        JMP         iC_RIi

iC_FRIi:SHUFPS      XMM7, XMM7, shift2LC            ;XMM7[0]=numComb[h]
        SHUFPS      XMM6, XMM6, shift2LC            ;XMM6[0]=c
        MOV         EDI, EDX                        ;EDI=deg
        SHL         EDI, 2                          ;EDI*4
        ADD         ECX, EDI                        ;ECX+deg*4, mi sposto di 1 riga su mComb
;        PEXTRD      EDI, XMM3, 00000011b            ;EDI=indice
        MOVSS       [EBX], XMM0                     ;xast[indice]=accumulatore
;        INC         EDI
        ADD         EBX, 4
;        PINSRD      XMM3, EDI, 00000011b            ;XMM3[0]=indice++
        ADDSS       XMM6, XMM5                      ;c++
        JMP         iC_RIc

;        extern void init_padd(VECTOR pDouble, MATRIX pDouble1, int n, int mt);

global init_padd

    iP_ar   equ 8
    iP_m    equ 12
    iP_mt   equ 16
    iP_n    equ 20

init_padd:
        start
        MOV         EAX, [EBP+iP_ar]
        MOVUPS      XMM0, [EAX]
        MOV         EAX, [EBP+iP_m]
        MOV         ECX, [EBP+iP_n]
        MOV         EBX, [EBP+iP_mt]
        SHL         EBX, 2
        DEC         ECX
        DEC         ECX

iP_c:   MOVUPS      [EAX], XMM0
        ADD         EAX, EBX
        LOOP        iP_c

        stop