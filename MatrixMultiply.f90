!Very rudimentary Fox Matrix Multiplication program
!@author Samuel Shadwell
!
!KNOWN ISSUES: Only works with square matrices
!   Requires more than 1 MPI process


PROGRAM MatrixMultiply
implicit none 
INCLUDE 'mpif.h'

!Variable declarations
INTEGER                                 :: ierr, commSize, rank                 !MPI variables
REAL, ALLOCATABLE, DIMENSION(:,:)       :: arayA, arayB, arayC                  !Matrices
INTEGER                                 :: arayLen, nodeRow, nodeCol
REAL                                    :: nodeAVal, nodeBVal, workA, workB
REAL                                    :: workC

!Initialize MPI Environment
CALL MPI_INIT(ierr)
CALL MPI_COMM_SIZE (MPI_COMM_WORLD, commSize, ierr)
CALL MPI_COMM_RANK (MPI_COMM_WORLD, rank, ierr)

!Allocates arrays
IF (rank .EQ. 0) THEN

    arayLen = 1
    ALLOCATE(arayA(arayLen, arayLen))
    ALLOCATE(arayB(arayLen, arayLen))
    ALLOCATE(arayC(arayLen, arayLen))

ELSE

    ALLOCATE(arayA(1,1))
    ALLOCATE(arayB(1,1))
    ALLOCATE(arayC(1,1))

ENDIF

!Sets and distributes array values to various nodes
CALL setUpArrays(ierr, commSize, rank, arayA, arayB, arayC, nodeAVal, nodeBVal, workA, workB, workC, arayLen)

!TODO: Write this code
!Passes values
!CALL passValues(?)

!Performs calculations
!CALL performCalculations(arayA, arayB, arayC, rank)

!Parallel Code Ends
CALL MPI_FINALIZE (ierr)

!Final Serial Code


END PROGRAM MatrixMultiply



!============== SUBROUTINE SetUpArrays ===============
!   Puts entire array on head node (0), and distributes single elements to 
!   all other nodes. Also sets size of arayA, B, and C to 1 for all
!   MPI processes other than the head node.
SUBROUTINE setUpArrays(ierr, commSize, rank, arayA, arayB, arayC, nodeAVal, nodeBVal, workA, workB, workC, arayLen)
implicit none 
include 'mpif.h'
INTEGER                                             :: counter,ierr, commSize, rank, arayLen, mpiStatus(MPI_STATUS_SIZE)
REAL, DIMENSION(arayLen,arayLen)                    :: arayA, arayB, arayC
REAL                                                :: nodeAVal, nodeBVal
REAL, INTENT(OUT)                                   :: workA, workB, workC
INTEGER                                             :: i, j, n, tag, nodeRow, nodeCol
REAL aOut, bOut

!Currently sets values for arrays, I need to figure out how to read them in

IF (rank .EQ. 0) THEN       ! Begin head node

    !Populates arrays A and B
    !TODO: replace this with an actual parser
    DO i = 1, arayLen
        DO j = 1, arayLen
            arayA(i,j) = i+j
            arayB(i,j) = i-j
            arayC(i,j) = 0  !Initialize all C elements as 0
        ENDDO
    ENDDO

    !Sends all nodes their values
    destinationRank = 0
    n = arayLen**2          !Number of data values which must be sent
    DO i = 1, arayLen
        DO j = 1, arayLen
            aOut = arayA(i,j)
            bOut = arayB(i,j)
            CALL MPI_SEND(aOut, 1, MPI_REAL, destinationRank, 1, MPI_COMM_WORLD, ierr)
           ! CALL MPI_SEND(bOut, 1, MPI_REAL, destinationRank, 2, MPI_COMM_WORLD, ierr)
            destinationRank = destinationRank + 1
        ENDDO
    ENDDO

ENDIF                       ! End head node


!Gets value for arrays A and B sent from head node


CALL MPI_RECV(nodeAVal, 1, MPI_REAL, 0, 1, MPI_COMM_WORLD, mpiStatus, ierr)
!CALL MPI_RECV(nodeBVal, 1, MPI_REAL, 0, 2, MPI_COMM_WORLD, mpiStatus, ierr)

END SUBROUTINE setUpArrays



SUBROUTINE performCalculations (arayA, arayB, arayC, rank)

REAL, ALLOCATABLE, DIMENSION(:,:)       :: arayA, arayB, arayC
INTEGER                                 :: rank

IF(rank .NE. 0) THEN
    
    !!arayC = arayC + (arayB*arayA)

ENDIF

END SUBROUTINE performCalculations
