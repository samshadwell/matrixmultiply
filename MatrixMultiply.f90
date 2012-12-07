!Very  rudimentary Fox Matrix Multiplication program
!@author Samuel Shadwell
!
!KNOWN ISSUES: Only works with 3x3 matrices
!   Limited modularity


PROGRAM MatrixMultiply
IMPLICIT NONE 
INCLUDE 'mpif.h'

!Variable declarations
INTEGER                                 :: ierr, commSize, worldRank
INTEGER                                 :: r1Comm, r2Comm, r3Comm, c1Comm, c2Comm, c3Comm !MPI variables
REAL, ALLOCATABLE, DIMENSION(:,:)       :: arayA, arayB, arayC                            !Matrices
REAL                                    :: nodeAVal, nodeBVal, workA, workB, workC
INTEGER, PARAMETER                      :: arayLen = 3 !Size of matrices to be multiplied (arayLen x arayLen)

!Initialize MPI Environment
CALL MPI_INIT(ierr)
CALL MPI_COMM_SIZE (MPI_COMM_WORLD, commSize, ierr)
CALL MPI_COMM_RANK (MPI_COMM_WORLD, worldRank, ierr)

!Allocates arrays
IF (worldRank .EQ. 0) THEN

    ALLOCATE(arayA(arayLen, arayLen))
    ALLOCATE(arayB(arayLen, arayLen))
    ALLOCATE(arayC(arayLen, arayLen))

!ELSE
!
!    ALLOCATE(arayA(1,1))
!    ALLOCATE(arayB(1,1))
!    ALLOCATE(arayC(1,1))

ENDIF

!Sets and distributes array values to various nodes
CALL setUpArrays(ierr, commSize, worldRank, arayA, arayB, arayC, nodeAVal, nodeBVal, arayLen)

!Creates communicators and groups for matrix rows and columns
CALL establishGroups(ierr, arayLen, r1Comm, r2Comm, r3Comm, c1Comm, c2Comm, c3Comm)

!Performs calculations
CALL performCalculations(ierr, workA, workB, workC, arayLen, r1Comm, r2Comm, r3Comm, c1Comm, c2Comm, c3Comm, nodeAVal, nodeBVal)

!Parallel Code Ends
CALL MPI_FINALIZE (ierr)

!Final Serial Code


END PROGRAM MatrixMultiply


!============== SUBROUTINE establishGroups ==========
!   Creates MPI groups to allow broadcasting of 
!	matrix values across the appropriate rows
!	and columns
!
!   NOTE: THIS IS NOT MODULAR DUE TO LIMITATIONS OF
!   MYSELF AND FORTRAN
!====================================================
SUBROUTINE establishGroups(ierr, arayLen, r1Comm, r2Comm, r3Comm, c1Comm, c2Comm, c3Comm)
IMPLICIT NONE
INCLUDE 'mpif.h'

INTEGER, INTENT(OUT)                               :: r1Comm, r2Comm, r3Comm   !Row communicators
INTEGER, INTENT(OUT)                               :: c1Comm, c2Comm, c3Comm   !Column communicators

INTEGER                                            :: r1Grp, r2Grp, r3Grp      !Row groups
INTEGER                                            :: c1Grp, c2Grp, c3Grp      !Column groups

INTEGER, INTENT(IN)                                :: arayLen
INTEGER, DIMENSION(3)                              :: r1Ranks, r2Ranks, r3Ranks
INTEGER, DIMENSION(3)                              :: c1Ranks, c2Ranks, c3Ranks
INTEGER                                            :: wGroup, i, ierr, testSize

CALL MPI_Comm_group(MPI_COMM_WORLD, wGroup)

!Establishes which worldRanks belong to each row and column, placing the
! indeces into *Ranks arrays
DO i = 0, arayLen-1

    r1Ranks(i+1) = i
    r2Ranks(i+1) = i + arayLen
    r3Ranks(i+1) = i + 2*arayLen

    c1Ranks(i+1) = i*arayLen
    c2Ranks(i+1) = i*arayLen + 1
    c3Ranks(i+1) = i*arayLen + 2

ENDDO

IF(ierr .NE. MPI_SUCCESS) THEN
    PRINT*, "error before group creation"
    STOP
ENDIF

!Creates 6 (arayLen*2) groups from the given worldRanks
CALL MPI_Group_incl(wGroup, arayLen, r1Ranks, r1Grp, ierr)
CALL MPI_Group_incl(wGroup, arayLen, r2Ranks, r2Grp, ierr)
CALL MPI_Group_incl(wGroup, arayLen, r3Ranks, r3Grp, ierr)
CALL MPI_Group_incl(wGroup, arayLen, c1Ranks, c1Grp, ierr)
CALL MPI_Group_incl(wGroup, arayLen, c2Ranks, c2Grp, ierr)
CALL MPI_Group_incl(wGroup, arayLen, c3Ranks, c3Grp, ierr) 
IF (ierr .NE. MPI_SUCCESS ) THEN
    PRINT*, "group creation error"
    STOP
ENDIF

!Creates communicators from the given groups
CALL MPI_Comm_create(MPI_COMM_WORLD, r1Grp, r1Comm , ierr)
CALL MPI_Comm_create(MPI_COMM_WORLD, r2Grp, r2Comm , ierr)
CALL MPI_Comm_create(MPI_COMM_WORLD, r3Grp, r3Comm , ierr)
CALL MPI_Comm_create(MPI_COMM_WORLD, c1Grp, c1Comm , ierr)
CALL MPI_Comm_create(MPI_COMM_WORLD, c2Grp, c2Comm , ierr)
CALL MPI_Comm_create(MPI_COMM_WORLD, c3Grp, c3Comm , ierr)
IF(ierr .NE. MPI_SUCCESS) THEN
    PRINT*, "communicator creation error"
    STOP
ENDIF


END SUBROUTINE establishGroups


!============== SUBROUTINE SetUpArrays ===============
!   Puts entire array on head node (0), and distributes single elements to 
!   all other nodes. Also sets size of arayA, B, and C to 1 for all
!   MPI processes other than the head node.
!====================================================
SUBROUTINE setUpArrays(ierr, commSize, worldRank, arayA, arayB, arayC, nodeAVal, nodeBVal, arayLen)
IMPLICIT NONE 
INCLUDE 'mpif.h'

INTEGER                                             :: counter, ierr, commSize, worldRank, arayLen, mpiStatus(MPI_STATUS_SIZE)
REAL, DIMENSION(arayLen,arayLen)                    :: arayA, arayB, arayC
REAL                                                :: nodeAVal, nodeBVal
INTEGER                                             :: i, j, tag, nodeRow, nodeCol, destinationRank
REAL aOut, bOut

!Currently sets values for arrays, I need to figure out how to read them in

IF (worldRank .EQ. 0) THEN       ! Begin head node

    !Populates arrays A and B
    DO i = 1, arayLen
        DO j = 1, arayLen
            arayA(i,j) = i+j
            arayB(i,j) = i-j
            arayC(i,j) = 0  !Initialize all C elements as 0
        ENDDO
    ENDDO

    !Sends all nodes their values
    destinationRank = 0
    DO i = 1, arayLen
        DO j = 1, arayLen
            IF (destinationRank .NE. 0) THEN
                aOut = arayA(i,j)
                bOut = arayB(i,j)
                CALL MPI_SEND(aOut, 1, MPI_REAL, destinationRank, 1, MPI_COMM_WORLD, ierr)
                CALL MPI_SEND(bOut, 1, MPI_REAL, destinationRank, 2, MPI_COMM_WORLD, ierr)
            ELSE
                nodeAVal = arayA(i,j)
                nodeBVal = arayB(i,j)
            ENDIF            
            destinationRank = destinationRank + 1
        ENDDO
    ENDDO

ENDIF                       ! End head node


!Gets value for arrays A and B sent from head node
IF (worldRank .NE. 0) THEN
    CALL MPI_RECV(nodeAVal, 1, MPI_REAL, 0, 1, MPI_COMM_WORLD, mpiStatus, ierr)
    CALL MPI_RECV(nodeBVal, 1, MPI_REAL, 0, 2, MPI_COMM_WORLD, mpiStatus, ierr)
ENDIF

END SUBROUTINE setUpArrays


!========================SUBROUTINE performCalculations==========
!   Subroutine which controls the broadcasting of values across
!   nodes as well as the calculating of C value
!================================================================
SUBROUTINE performCalculations(ierr, workA, workB, workC, arayLen, r1Comm, r2Comm, r3Comm, c1Comm, c2Comm, c3Comm, nodeAVal, nodeBVal)
INTEGER                                            :: arayLen, i, r1Comm, r2Comm, r3Comm, c1Comm, c2Comm, c3Comm, ierr
REAL                                               :: workA, workB, workC, nodeAVal, nodeBVal

workC = 0.0d0

DO i = 0, arayLen-1

    workA = nodeAVal
    workB = nodeBVal
    CALL passValues(i, workA, workB, workC, r1Comm, r2Comm, r3Comm, c1Comm, c2Comm, c3Comm, ierr)
    CALL calculate(workA, workB, workC)

ENDDO

OPEN(UNIT = 11, FILE = "finalMatrix.txt")
WRITE(11,*) "NODE: ", worldRank, "   VALUE: ", workC
CLOSE(11)

END SUBROUTINE performCalculations


!========================SUBROUTINE calculate====================
!   Relatively simple subroutine which performs the value addition
!   and multiplication (C = C + AxB)
!================================================================
SUBROUTINE calculate(workA, workB, workC)
REAL                                                :: workA, workB, workC

workC = workC + (workB*workA)

END SUBROUTINE calculate


!========================SUBROUTINE passValues===================
!   Subroutine which passes the necessary values to the arrays
!================================================================
SUBROUTINE passValues(i, workA, workB, workC, r1Comm, r2Comm, r3Comm, c1Comm, c2Comm, c3Comm)
INCLUDE 'mpif.h'
INTEGER                                            :: rank, i, r1Comm, r2Comm, r3Comm, c1Comm, c2Comm, c3Comm, ierr
REAL                                               :: workA, workB, workC

CALL MPI_COMM_RANK (MPI_COMM_WORLD, rank, ierr)

IF ((rank .LE. 2) .AND. (rank .GE. 0)) THEN
    CALL MPI_Bcast(workA, 1, MPI_REAL, i, r1Comm, ierr)
ELSEIF((rank .LE. 5) .AND. (rank .GE. 3)) THEN
    CALL MPI_Bcast(workA, 1, MPI_REAL, i, r2Comm, ierr)
ELSEIF((rank .LE. 8) .AND. (rank .GE. 6)) THEN
    CALL MPI_Bcast(workA, 1, MPI_REAL, i, r3Comm, ierr)
ENDIF

IF    ((rank .EQ. 0) .OR. (rank .EQ. 3) .OR. (rank .EQ. 6)) THEN
    CALL MPI_Bcast(workB, 1, MPI_REAL, i, c1Comm, ierr)
ELSEIF((rank .EQ. 1) .OR. (rank .EQ. 4) .OR. (rank .EQ. 7)) THEN
    CALL MPI_Bcast(workB, 1, MPI_REAL, i, c2Comm, ierr)
ELSEIF((rank .EQ. 2) .OR. (rank .EQ. 5) .OR. (rank .EQ. 8)) THEN
    CALL MPI_Bcast(workB, 1, MPI_REAL, i, c3Comm, ierr)
ENDIF

END SUBROUTINE passValues
