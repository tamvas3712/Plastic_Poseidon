    integer :: pii,o2o,o3c,o4n,n1p,n3n,n4n,n5s,r6c,r6n,r6p,r6s,r1c,r1n, &
    r1p,p1c,p1n,p1p,p1s,p2c,p2n,p2p,p3c,p3n,p3p,p4c,p4n,p4p, &
    z4c,z5c,z5n,z5p,z6c,z6n,z6p,b1c,b1n,b1p,aii,a1c,a1n,a1p, &
    a1s,h1c,h2c,q1c,q1n,q1p,q6c,q6n,q6p,q6s,q7c,q7n,q7p,q17c, &
    q17n,q17p,g2o,g3c,g4n,k5s,k15s,k25s,k6e,k26e,k4n,k14n, &
    k24n,k1p,k11p,k21p,k3n,k13n,k23n,d1m,d2m,xd2m,d3m,d4m, &
    d5m,d6m,d7m,d8m,d9m,y2c,y3c,y4c,talk

    parameter &
!      Pelagic state variables
    (pii=1, o2o=2, o3c=3, o4n=4, n1p=5, n3n=6, n4n=7, &
    n5s=8, r6c=9, r6n=10, r6p=11, r6s=12, r1c=13, &
    r1n=14, r1p=15, p1c=16, p1n=17, p1p=18, p1s=19, &
    p2c=20, p2n=21, p2p=22, p3c=23, p3n=24, p3p=25, &
    p4c=26, p4n=27, p4p=28, z4c=29, z5c=30, z5n=31, &
    z5p=32, z6c=33, z6n=34, z6p=35, b1c=36, b1n=37, &
    b1p=38,talk=39, &
!      Benthic state variables
    aii=40, a1c=41, a1n=42, a1p=43, a1s=44, &
    h1c=45, h2c=46, q1c=47, q1n=48, q1p=49, q6c=50, &
    q6n=51, q6p=52, q6s=53, q7c=54, q7n=55, q7p=56, &
    q17c=57, q17n=58, q17p=59, g2o=60, g3c=61, g4n=62, &
    k5s=63, k15s=64, k25s=65, k6e=66, k26e=67, k4n=68, &
    k14n=69, k24n=70, k1p=71, k11p=72, k21p=73, k3n=74, &
    k13n=75, k23n=76, d1m=77, d2m=78, xd2m=79, d3m=80, &
    d4m=81, d5m=82, d6m=83, d7m=84, d8m=85, d9m=86, y2c=87, &
    y3c=88, y4c=89)



    INTEGER :: P_STATE, I_STATE, B_STATE,N_COMPg, n_upper$g
    INTEGER :: N_COMP, n_upper$,N_COMPben
! How often to solve the ecology in relation with the physics
    integer :: ecosteps
    parameter (ecosteps = 1)
                        	
! for the assimilation filter
    INTEGER ::   N_STATE,I_COMPb

! n_upper$ = N_COMP/kbm1
    PARAMETER (I_STATE=89, N_COMPg=610392, n_upper$g = 25433)
    PARAMETER(P_STATE=39,B_STATE=50)

    PARAMETER (N_COMPben=n_upper$g)
!      PARAMETER(P_STATE=6,B_STATE=3)

    PARAMETER (N_COMP=N_COMPg,n_upper$=n_upper$g)
    PARAMETER (I_COMPb=N_COMP-n_upper$+1)
!      integer N_COMP,n_upper$
!      COMMON /BOECOSIZE/ N_COMP,n_upper$
