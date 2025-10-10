module m_spectral

!------------------------------------------------------------------------------------!
!! Purpose:                                                                         -!
!  Contains the all necessary information regarding spectral elements including:    -!
!  - Lobatto polynomials                                                            -!
!  - Legendre polynomials                                                           -!
!  - Lobatto quadrature points and integration weights                              -!
!                                                                                   -!
!! Record of revisions:                                                             -!
!   Date       Programmer     Description of change                                 -!
!   ====       ==========     =====================                                 -!
! 04/09/24      C. Jones         Original code                                      -!
!------------------------------------------------------------------------------------!

use m_constants
use m_utilib
use m_quadrature
implicit none

type :: t_spectral
    real(dp), dimension(:), allocatable :: Xi, W
end type t_spectral

contains 

subroutine Spectral1DPositions(Spectral, PolyOrder)
    type(t_spectral), intent(inout)     :: Spectral
    integer, intent(in)                 :: PolyOrder

    allocate(Spectral%Xi(PolyOrder+1))

    Spectral%Xi(1) = -1.0_dp
    Spectral%Xi(PolyOrder+1) = 1.0_dp
    
    select case (PolyOrder)
    case(2)
        Spectral%Xi(2) = 0.0_dp
    case(3)
        Spectral%Xi(2) = -1.0_dp / sqrt(5.0_dp)
        Spectral%Xi(3) = -Spectral%Xi(2)
    case(4)
        Spectral%Xi(2) = -sqrt(3.0_dp/7.0_dp)
        Spectral%Xi(3) = 0.0_dp
        Spectral%Xi(4) = -Spectral%Xi(2)
    case(5)
        Spectral%Xi(2) = -0.765055323929465_dp
        Spectral%Xi(3) = -0.285231516480645_dp
        Spectral%Xi(4) = -Spectral%Xi(3)
        Spectral%Xi(5) = -Spectral%Xi(2)
    case(6)
        Spectral%Xi(2) = -0.830223896278567_dp
        Spectral%Xi(3) = -0.468848793470714_dp
        Spectral%Xi(4) = 0.0_dp
        Spectral%Xi(5) = -Spectral%Xi(3)
        Spectral%Xi(6) = -Spectral%Xi(2)
    case(7)
        Spectral%Xi(2) = -0.871740148509606_dp
        Spectral%Xi(3) = -0.591700181433142_dp
        Spectral%Xi(4) = -0.209299217902478_dp
        Spectral%Xi(5) = -Spectral%Xi(4)
        Spectral%Xi(6) = -Spectral%Xi(3)
        Spectral%Xi(7) = -Spectral%Xi(2)
    case(8)
        Spectral%Xi(2) = -0.899757995411460_dp
        Spectral%Xi(3) = -0.677186279510737_dp
        Spectral%Xi(4) = -0.363117463826178_dp
        Spectral%Xi(5) = 0.0_dp
        Spectral%Xi(6) = -Spectral%Xi(4)
        Spectral%Xi(7) = -Spectral%Xi(3)
        Spectral%Xi(8) = -Spectral%Xi(2)
    case(9)
        Spectral%Xi(2) = -0.919533908166459_dp
        Spectral%Xi(3) = -0.738773865105505_dp
        Spectral%Xi(4) = -0.477924949810444_dp
        Spectral%Xi(5) = -0.165278957666387_dp
        Spectral%Xi(6) = -Spectral%Xi(5)
        Spectral%Xi(7) = -Spectral%Xi(4)
        Spectral%Xi(8) = -Spectral%Xi(3)
        Spectral%Xi(9) = -Spectral%Xi(2)
    case(10)
        Spectral%Xi(2) = -0.934001430408059_dp
        Spectral%Xi(3) = -0.784483473663144_dp
        Spectral%Xi(4) = -0.565235326996205_dp
        Spectral%Xi(5) = -0.295758135586939_dp
        Spectral%Xi(6) = 0.0_dp
        Spectral%Xi(7) = -Spectral%Xi(5)
        Spectral%Xi(8) = -Spectral%Xi(4)
        Spectral%Xi(9) = -Spectral%Xi(3)
        Spectral%Xi(10) = -Spectral%Xi(2)

    case(11)
        Spectral%Xi(2) = -0.9448992722296681_dp
        Spectral%Xi(3) = -0.8192793216440067_dp
        Spectral%Xi(4) = -0.6328761530318606_dp
        Spectral%Xi(5) = -0.3995309409653489_dp
        Spectral%Xi(6) = -0.1365529328549276_dp
        Spectral%Xi(7) = -Spectral%Xi(6)
        Spectral%Xi(8) = -Spectral%Xi(5)
        Spectral%Xi(9) = -Spectral%Xi(4)
        Spectral%Xi(10) = -Spectral%Xi(3)
        Spectral%Xi(11) = -Spectral%Xi(2)

    case(12)
        Spectral%Xi(2) = -0.9533098466421639_dp
        Spectral%Xi(3) = -0.8463475646518723_dp
        Spectral%Xi(4) = -0.6861884690817574_dp
        Spectral%Xi(5) = -0.4829098210913362_dp
        Spectral%Xi(6) = -0.249286930106240_dp
        Spectral%Xi(7) = 0.0_dp
        Spectral%Xi(8) = -Spectral%Xi(6)
        Spectral%Xi(9) = -Spectral%Xi(5)
        Spectral%Xi(10) = -Spectral%Xi(4)
        Spectral%Xi(11) = -Spectral%Xi(3)
        Spectral%Xi(12) = -Spectral%Xi(2)

    case(13)
        Spectral%Xi(2) = -0.959935045267261_dp
        Spectral%Xi(3) = -0.867801053830347_dp
        Spectral%Xi(4) = -0.728868599091326_dp
        Spectral%Xi(5) = -0.550639402928647_dp
        Spectral%Xi(6) = -0.342724013342712_dp
        Spectral%Xi(7) = -0.116331868883703_dp
        Spectral%Xi(8) = -Spectral%Xi(7)
        Spectral%Xi(9) = -Spectral%Xi(6)
        Spectral%Xi(10) = -Spectral%Xi(5)
        Spectral%Xi(11) = -Spectral%Xi(4)
        Spectral%Xi(12) = -Spectral%Xi(3)
        Spectral%Xi(13) = -Spectral%Xi(2)

    case(14)
        Spectral%Xi(2) = -0.965245926503839_dp
        Spectral%Xi(3) = -0.885082044222976_dp
        Spectral%Xi(4) = -0.763519689951815_dp
        Spectral%Xi(5) = -0.606253205469845_dp
        Spectral%Xi(6) = -0.420638054713672_dp
        Spectral%Xi(7) = -0.215353955363794_dp
        Spectral%Xi(8) = 0.0_dp
        Spectral%Xi(9) = -Spectral%Xi(7)
        Spectral%Xi(10) = -Spectral%Xi(6)
        Spectral%Xi(11) = -Spectral%Xi(5)
        Spectral%Xi(12) = -Spectral%Xi(4)
        Spectral%Xi(13) = -Spectral%Xi(3)
        Spectral%Xi(14) = -Spectral%Xi(2)

    case(15)
        Spectral%Xi(2) = -0.969568046270218_dp
        Spectral%Xi(3) = -0.899200533093472_dp
        Spectral%Xi(4) = -0.7920082918618151_dp
        Spectral%Xi(5) = -0.652388702882493_dp
        Spectral%Xi(6) = -0.486059421887137_dp
        Spectral%Xi(7) = -0.2998304689007632_dp
        Spectral%Xi(8) = -0.1013262735219491_dp
        Spectral%Xi(9) = -Spectral%Xi(8)
        Spectral%Xi(10) = -Spectral%Xi(7)
        Spectral%Xi(11) = -Spectral%Xi(6)
        Spectral%Xi(12) = -Spectral%Xi(5)
        Spectral%Xi(13) = -Spectral%Xi(4)
        Spectral%Xi(14) = -Spectral%Xi(3)
        Spectral%Xi(15) = -Spectral%Xi(2)

    case(16)
        Spectral%Xi(2) = -0.9731321766314183_dp
        Spectral%Xi(3) = -0.910879995915574_dp
        Spectral%Xi(4) = -0.8156962512217703_dp
        Spectral%Xi(5) = -0.6910289806276847_dp
        Spectral%Xi(6) = -0.541385399330102_dp
        Spectral%Xi(7) = -0.3721744335654772_dp
        Spectral%Xi(8) = -0.189511973518317_dp
        Spectral%Xi(9) = 0.0_dp
        Spectral%Xi(10) = -Spectral%Xi(8)
        Spectral%Xi(11) = -Spectral%Xi(7)
        Spectral%Xi(12) = -Spectral%Xi(6)
        Spectral%Xi(13) = -Spectral%Xi(5)
        Spectral%Xi(14) = -Spectral%Xi(4)
        Spectral%Xi(15) = -Spectral%Xi(3)
        Spectral%Xi(16) = -Spectral%Xi(2)

    case(17)
        Spectral%Xi(2) = -0.976105557412198_dp
        Spectral%Xi(3) = -0.920649185347533_dp
        Spectral%Xi(4) = -0.835593535218090_dp
        Spectral%Xi(5) = -0.723679329283243_dp
        Spectral%Xi(6) = -0.588504834318661_dp
        Spectral%Xi(7) = -0.434415036912123_dp
        Spectral%Xi(8) = -0.2663626528782805_dp
        Spectral%Xi(9) = -0.089749093484652_dp
        Spectral%Xi(10) = -Spectral%Xi(9)
        Spectral%Xi(11) = -Spectral%Xi(8)
        Spectral%Xi(12) = -Spectral%Xi(7)
        Spectral%Xi(13) = -Spectral%Xi(6)
        Spectral%Xi(14) = -Spectral%Xi(5)
        Spectral%Xi(15) = -Spectral%Xi(4)
        Spectral%Xi(16) = -Spectral%Xi(3)
        Spectral%Xi(17) = -Spectral%Xi(2)

    case(18)
        Spectral%Xi(2) = -0.978611766222080_dp
        Spectral%Xi(3) = -0.928901528152586_dp
        Spectral%Xi(4) = -0.852460577796646_dp
        Spectral%Xi(5) = -0.751494202552613_dp
        Spectral%Xi(6) = -0.628908137265221_dp
        Spectral%Xi(7) = -0.488229285680714_dp
        Spectral%Xi(8) = -0.333504847824499_dp
        Spectral%Xi(9) = -0.169186023409282_dp
        Spectral%Xi(10) = 0.0_dp
        Spectral%Xi(11) = -Spectral%Xi(9)
        Spectral%Xi(12) = -Spectral%Xi(8)
        Spectral%Xi(13) = -Spectral%Xi(7)
        Spectral%Xi(14) = -Spectral%Xi(6)
        Spectral%Xi(15) = -Spectral%Xi(5)
        Spectral%Xi(16) = -Spectral%Xi(4)
        Spectral%Xi(17) = -Spectral%Xi(3)
        Spectral%Xi(18) = -Spectral%Xi(2)

    case(19)
        Spectral%Xi(2) = -0.980743704893914_dp
        Spectral%Xi(3) = -0.935934498812665_dp
        Spectral%Xi(4) = -0.866877978089950_dp
        Spectral%Xi(5) = -0.775368260952056_dp
        Spectral%Xi(6) = -0.663776402290311_dp
        Spectral%Xi(7) = -0.534992864031886_dp
        Spectral%Xi(8) = -0.392353183713909_dp
        Spectral%Xi(9) = -0.239551705922986_dp
        Spectral%Xi(10) = -0.080545937238822_dp
        Spectral%Xi(11) = -Spectral%Xi(10)
        Spectral%Xi(12) = -Spectral%Xi(9)
        Spectral%Xi(13) = -Spectral%Xi(8)
        Spectral%Xi(14) = -Spectral%Xi(7)
        Spectral%Xi(15) = -Spectral%Xi(6)
        Spectral%Xi(16) = -Spectral%Xi(5)
        Spectral%Xi(17) = -Spectral%Xi(4)
        Spectral%Xi(18) = -Spectral%Xi(3)
        Spectral%Xi(19) = -Spectral%Xi(2)

    case(20)
        Spectral%Xi(2) = -0.982572296604548_dp
        Spectral%Xi(3) = -0.941976296959746_dp
        Spectral%Xi(4) = -0.879294755323591_dp
        Spectral%Xi(5) = -0.796001926077712_dp
        Spectral%Xi(6) = -0.694051026062223_dp
        Spectral%Xi(7) = -0.575831960261831_dp
        Spectral%Xi(8) = -0.444115783279002_dp
        Spectral%Xi(9) = -0.301989856508765_dp
        Spectral%Xi(10) = -0.152785515802186_dp
        Spectral%Xi(11) = 0.0_dp
        Spectral%Xi(12) = -Spectral%Xi(10)
        Spectral%Xi(13) = -Spectral%Xi(9)
        Spectral%Xi(14) = -Spectral%Xi(8)
        Spectral%Xi(15) = -Spectral%Xi(7)
        Spectral%Xi(16) = -Spectral%Xi(6)
        Spectral%Xi(17) = -Spectral%Xi(5)
        Spectral%Xi(18) = -Spectral%Xi(4)
        Spectral%Xi(19) = -Spectral%Xi(3)
        Spectral%Xi(20) = -Spectral%Xi(2)

    end select
    

end subroutine Spectral1DPositions

subroutine SE_Pos(Xi, IntegOrder)
    real(dp), dimension(:), allocatable :: Xi
    integer, intent(in)                 :: IntegOrder   

    allocate(Xi(IntegOrder+1))

    Xi(1) = -1.0_dp
    Xi(IntegOrder+1) = 1.0_dp
    
    select case (IntegOrder)
    case(2)
        Xi(2) = 0.0_dp
    case(3)
        Xi(2) = -1.0_dp / sqrt(5.0_dp)
        Xi(3) = -Xi(2)
    case(4)
        Xi(2) = -sqrt(3.0_dp/7.0_dp)
        Xi(3) = 0.0_dp
        Xi(4) = -Xi(2)
    case(5)
        Xi(2) = -0.765055323929465_dp
        Xi(3) = -0.285231516480645_dp
        Xi(4) = -Xi(3)
        Xi(5) = -Xi(2)
    case(6)
        Xi(2) = -0.830223896278567_dp
        Xi(3) = -0.468848793470714_dp
        Xi(4) = 0.0_dp
        Xi(5) = -Xi(3)
        Xi(6) = -Xi(2)
    case(7)
        Xi(2) = -0.871740148509606_dp
        Xi(3) = -0.591700181433142_dp
        Xi(4) = -0.209299217902478_dp
        Xi(5) = -Xi(4)
        Xi(6) = -Xi(3)
        Xi(7) = -Xi(2)
    case(8)
        Xi(2) = -0.899757995411460_dp
        Xi(3) = -0.677186279510737_dp
        Xi(4) = -0.363117463826178_dp
        Xi(5) = 0.0_dp
        Xi(6) = -Xi(4)
        Xi(7) = -Xi(3)
        Xi(8) = -Xi(2)
    case(9)
        Xi(2) = -0.919533908166459_dp
        Xi(3) = -0.738773865105505_dp
        Xi(4) = -0.477924949810444_dp
        Xi(5) = -0.165278957666387_dp
        Xi(6) = -Xi(5)
        Xi(7) = -Xi(4)
        Xi(8) = -Xi(3)
        Xi(9) = -Xi(2)
    case(10)
        Xi(2) = -0.934001430408059_dp
        Xi(3) = -0.784483473663144_dp
        Xi(4) = -0.565235326996205_dp
        Xi(5) = -0.295758135586939_dp
        Xi(6) = 0.0_dp
        Xi(7) = -Xi(5)
        Xi(8) = -Xi(4)
        Xi(9) = -Xi(3)
        Xi(10) = -Xi(2)

    case(11)
        Xi(2) = -0.9448992722296681_dp
        Xi(3) = -0.8192793216440067_dp
        Xi(4) = -0.6328761530318606_dp
        Xi(5) = -0.3995309409653489_dp
        Xi(6) = -0.1365529328549276_dp
        Xi(7) = -Xi(6)
        Xi(8) = -Xi(5)
        Xi(9) = -Xi(4)
        Xi(10) = -Xi(3)
        Xi(11) = -Xi(2)

    case(12)
        Xi(2) = -0.9533098466421639_dp
        Xi(3) = -0.8463475646518723_dp
        Xi(4) = -0.6861884690817574_dp
        Xi(5) = -0.4829098210913362_dp
        Xi(6) = -0.249286930106240_dp
        Xi(7) = 0.0_dp
        Xi(8) = -Xi(6)
        Xi(9) = -Xi(5)
        Xi(10) = -Xi(4)
        Xi(11) = -Xi(3)
        Xi(12) = -Xi(2)

    case(13)
        Xi(2) = -0.959935045267261_dp
        Xi(3) = -0.867801053830347_dp
        Xi(4) = -0.728868599091326_dp
        Xi(5) = -0.550639402928647_dp
        Xi(6) = -0.342724013342712_dp
        Xi(7) = -0.116331868883703_dp
        Xi(8) = -Xi(7)
        Xi(9) = -Xi(6)
        Xi(10) = -Xi(5)
        Xi(11) = -Xi(4)
        Xi(12) = -Xi(3)
        Xi(13) = -Xi(2)

    case(14)
        Xi(2) = -0.965245926503839_dp
        Xi(3) = -0.885082044222976_dp
        Xi(4) = -0.763519689951815_dp
        Xi(5) = -0.606253205469845_dp
        Xi(6) = -0.420638054713672_dp
        Xi(7) = -0.215353955363794_dp
        Xi(8) = 0.0_dp
        Xi(9) = -Xi(7)
        Xi(10) = -Xi(6)
        Xi(11) = -Xi(5)
        Xi(12) = -Xi(4)
        Xi(13) = -Xi(3)
        Xi(14) = -Xi(2)

    case(15)
        Xi(2) = -0.969568046270218_dp
        Xi(3) = -0.899200533093472_dp
        Xi(4) = -0.7920082918618151_dp
        Xi(5) = -0.652388702882493_dp
        Xi(6) = -0.486059421887137_dp
        Xi(7) = -0.2998304689007632_dp
        Xi(8) = -0.1013262735219491_dp
        Xi(9) = -Xi(8)
        Xi(10) = -Xi(7)
        Xi(11) = -Xi(6)
        Xi(12) = -Xi(5)
        Xi(13) = -Xi(4)
        Xi(14) = -Xi(3)
        Xi(15) = -Xi(2)

    case(16)
        Xi(2) = -0.9731321766314183_dp
        Xi(3) = -0.910879995915574_dp
        Xi(4) = -0.8156962512217703_dp
        Xi(5) = -0.6910289806276847_dp
        Xi(6) = -0.541385399330102_dp
        Xi(7) = -0.3721744335654772_dp
        Xi(8) = -0.189511973518317_dp
        Xi(9) = 0.0_dp
        Xi(10) = -Xi(8)
        Xi(11) = -Xi(7)
        Xi(12) = -Xi(6)
        Xi(13) = -Xi(5)
        Xi(14) = -Xi(4)
        Xi(15) = -Xi(3)
        Xi(16) = -Xi(2)

    case(17)
        Xi(2) = -0.976105557412198_dp
        Xi(3) = -0.920649185347533_dp
        Xi(4) = -0.835593535218090_dp
        Xi(5) = -0.723679329283243_dp
        Xi(6) = -0.588504834318661_dp
        Xi(7) = -0.434415036912123_dp
        Xi(8) = -0.2663626528782805_dp
        Xi(9) = -0.089749093484652_dp
        Xi(10) = -Xi(9)
        Xi(11) = -Xi(8)
        Xi(12) = -Xi(7)
        Xi(13) = -Xi(6)
        Xi(14) = -Xi(5)
        Xi(15) = -Xi(4)
        Xi(16) = -Xi(3)
        Xi(17) = -Xi(2)

    case(18)
        Xi(2) = -0.978611766222080_dp
        Xi(3) = -0.928901528152586_dp
        Xi(4) = -0.852460577796646_dp
        Xi(5) = -0.751494202552613_dp
        Xi(6) = -0.628908137265221_dp
        Xi(7) = -0.488229285680714_dp
        Xi(8) = -0.333504847824499_dp
        Xi(9) = -0.169186023409282_dp
        Xi(10) = 0.0_dp
        Xi(11) = -Xi(9)
        Xi(12) = -Xi(8)
        Xi(13) = -Xi(7)
        Xi(14) = -Xi(6)
        Xi(15) = -Xi(5)
        Xi(16) = -Xi(4)
        Xi(17) = -Xi(3)
        Xi(18) = -Xi(2)

    case(19)
        Xi(2) = -0.980743704893914_dp
        Xi(3) = -0.935934498812665_dp
        Xi(4) = -0.866877978089950_dp
        Xi(5) = -0.775368260952056_dp
        Xi(6) = -0.663776402290311_dp
        Xi(7) = -0.534992864031886_dp
        Xi(8) = -0.392353183713909_dp
        Xi(9) = -0.239551705922986_dp
        Xi(10) = -0.080545937238822_dp
        Xi(11) = -Xi(10)
        Xi(12) = -Xi(9)
        Xi(13) = -Xi(8)
        Xi(14) = -Xi(7)
        Xi(15) = -Xi(6)
        Xi(16) = -Xi(5)
        Xi(17) = -Xi(4)
        Xi(18) = -Xi(3)
        Xi(19) = -Xi(2)

    case(20)
        Xi(2) = -0.982572296604548_dp
        Xi(3) = -0.941976296959746_dp
        Xi(4) = -0.879294755323591_dp
        Xi(5) = -0.796001926077712_dp
        Xi(6) = -0.694051026062223_dp
        Xi(7) = -0.575831960261831_dp
        Xi(8) = -0.444115783279002_dp
        Xi(9) = -0.301989856508765_dp
        Xi(10) = -0.152785515802186_dp
        Xi(11) = 0.0_dp
        Xi(12) = -Xi(10)
        Xi(13) = -Xi(9)
        Xi(14) = -Xi(8)
        Xi(15) = -Xi(7)
        Xi(16) = -Xi(6)
        Xi(17) = -Xi(5)
        Xi(18) = -Xi(4)
        Xi(19) = -Xi(3)
        Xi(20) = -Xi(2)

    end select

end subroutine SE_Pos

real(dp) function Lobatto1D(PolyOrder, Xi) result (Lobatto)

    integer, intent(in)     :: PolyOrder
    real(dp), intent(in)    :: Xi

    lobatto = 0.0_dp
    select case (PolyOrder)
    case(0)
        Lobatto = 1.0
    case(1)
        Lobatto = 3.0_dp * Xi
    case(2)
        Lobatto = 1.5_dp * (5.0_dp * Xi**2.0_dp - 1.0_dp)
    case(3)
        Lobatto = 2.5_dp * (7.0_dp * Xi**2.0_dp - 3.0_dp) * Xi
    case(4)
        Lobatto = 1.875_dp * (21.0_dp * Xi**4.0_dp - 14.0_dp * Xi**2.0_dp + 1.0_dp)
    case(5)
        Lobatto = 0.125_dp * (693.0_dp * Xi**4.0_dp - 630.0_dp * Xi**2.0_dp + 105.0_dp) * Xi
    case(6)
        Lobatto = 0.0625_dp * (3003.0_dp * Xi**6.0_dp - 3465.0_dp * Xi**4.0_dp + 945.0_dp * Xi**2.0_dp - 35.0_dp)
    case(7)
        Lobatto = (9.0_dp / 16.0_dp) * Xi * (715.0_dp * Xi**6.0_dp - 1001.0_dp * Xi**4.0_dp + 385.0_dp * Xi**2.0_dp - 35.0_dp)
    case(8)
        Lobatto = (45.0_dp / 128.0_dp) * (2431.0_dp * Xi**8.0_dp - 4004.0_dp * Xi**6.0_dp + 2002.0_dp * Xi**4.0_dp - 308.0_dp * Xi**2.0_dp + 7.0_dp)
    case(9)
        Lobatto = (55.0_dp / 128.0_dp) * Xi * (4199.0_dp * Xi**8.0_dp - 7956.0_dp * Xi**6.0_dp + 4914.0_dp * Xi**4.0_dp - 1092.0_dp * Xi**2.0_dp + 63.0_dp)
    case(10)
        Lobatto = (33.0_dp / 256.0_dp) * (29393.0_dp * Xi**10.0_dp - 62985.0_dp * Xi**8.0_dp + 46410.0_dp * Xi**6.0_dp - 13650.0_dp * Xi**4.0_dp + 1365.0_dp * Xi**2.0_dp - 21.0_dp)
    end select

end function Lobatto1D

real(dp) function dLobatto_dx(PolyOrder, Xi) result(Val)

    integer, intent(in)     :: PolyOrder
    real(dp), intent(in)    :: Xi

    Val = 0.0_dp
    select case (PolyOrder)
    case(0)
        Val = 0.0_dp
    case(1)
        Val = 3.0_dp
    case(2)
        Val = 15.0_dp * Xi
    case(3)
        Val = (105.0_dp * Xi**2.0_dp - 15.0_dp) / 2.0_dp
    case(4)
        Val = (315.0_dp * Xi**3.0_dp - 105.0_dp * Xi) / 2.0_dp
    case(5)
        Val = (105.0_dp / 8.0_dp) * (1.0_dp - 18.0_dp * Xi**2.0_dp + 33.0_dp * Xi**4.0_dp) 
    case(6)
        Val = (63.0_dp / 8.0_dp) * Xi * (15.0_dp - 110.0_dp * Xi**2.0_dp + 143.0_dp * Xi**4.0_dp)
    case(7)
        Val = (315.0_dp / 16.0_dp) * (-1.0_dp + 33.0_dp * Xi**2.0_dp - 143.0_dp * Xi**4.0_dp + 143.0_dp * Xi**6.0_dp)
    case(8)
        Val = (495.0_dp / 16.0_dp) * Xi * (-7.0_dp + 91.0_dp * Xi**2.0_dp - 273.0_dp * Xi**4.0_dp + 221.0_dp * Xi**6.0_dp)
    case(9)
        Val = (495.0_dp / 128.0_dp) * (7.0_dp + 13.0_dp * Xi**2.0_dp * (-28.0_dp + 210.0_dp * Xi**2.0_dp - 476.0_dp * Xi**4.0_dp + 323.0_dp * Xi**6.0_dp))
    case(10)
        Val = (2145.0_dp / 128.0_dp) * Xi * (21.0_dp - 420.0_dp * Xi**2.0_dp + 2142.0_dp * Xi**4.0_dp - 3876.0_dp * Xi**6.0_dp + 2261.0_dp * Xi**8.0_dp)
    end select

end function dLobatto_dx

real(dp) function Legendre(PolyOrder, x) result (Val)

    integer, intent(in)     :: PolyOrder
    real(dp), intent(in)    :: x

    Val = 0.0_dp
    select case (PolyOrder)
    case(0)
        Val = 1.0_dp
    case(1)
        Val = x
    case(2)
        Val = 0.5_dp * (3.0_dp * x**2.0_dp - 1.0_dp)
    case(3)
        Val = 0.5_dp * (5.0_dp * x**3.0_dp - 3.0_dp * x)
    case(4)
        Val = 1.0_dp / 8.0_dp * (35.0_dp * x**4.0_dp - 30.0_dp * x**2.0_dp + 3.0_dp)
    case(5)
        Val = 1.0_dp / 8.0_dp * (63.0_dp * x**5.0_dp - 70.0_dp * x**3.0_dp + 15.0_dp * x)  
    case(6)
        Val = 1.0_dp / 16.0_dp * (231.0_dp * x**6.0_dp - 315.0_dp * x**4.0_dp + 105.0_dp * x**2.0_dp - 5.0_dp)
    case(7)
        Val = 1.0_dp / 16.0_dp * (429.0_dp * x**7.0_dp - 693.0_dp * x**5.0_dp + 315.0_dp * x**3.0_dp - 35.0_dp * x)
    case(8)
        Val = 1.0_dp / 128.0_dp * (6435.0_dp * x**8.0_dp - 12012.0_dp * x**6.0_dp + 6930.0_dp * x**4.0_dp - 1260.0_dp * x**2.0_dp + 35.0_dp)
    case(9)
        Val = 1.0_dp / 128.0_dp * (12155.0_dp * x**9.0_dp - 25740.0_dp * x**7.0_dp + 18018.0_dp * x**5.0_dp - 4620.0_dp * x**3.0_dp + 315.0_dp * x)
    case(10)
        Val = 1.0_dp / 256.0_dp * (46189.0_dp * x**10.0_dp - 109395.0_dp * x**8.0_dp + 90090.0_dp * x**6.0_dp - 30030.0_dp * x**4.0_dp + 3465.0_dp * x**2.0_dp - 63.0_dp)
    end select

end function Legendre

real(dp) function dLegendre_dx(PolyOrder, x) result(Val)

    integer, intent(in)     :: PolyOrder
    real(dp), intent(in)    :: x

    Val = 0.0_dp
    select case (PolyOrder)
    case(0)
        Val = 0.0_dp
    case(1)
        Val = 1.0_dp
    case(2)
        Val = 3.0_dp * x
    case(3)
        Val = (15.0_dp * x**2.0_dp - 3.0_dp) / 2.0_dp
    case(4)
        Val = (140.0_dp * x**3.0_dp - 60.0_dp * x) / 8.0_dp
    case(5)
        Val = (315.0_dp * x**4.0_dp - 210.0_dp * x**2.0_dp + 15.0_dp) / 8.0_dp
    case(6)
        Val = (1386.0_dp * x**5.0_dp - 1260.0_dp * x**3.0_dp + 210.0_dp * x) / 16.0_dp
    case(7)
        Val = (2574.0_dp * x**5.0_dp - 3465.0_dp * x**4.0_dp + 945.0_dp * x**2.0_dp - 35.0_dp) / 16.0_dp
    case(8)
        Val = (51480.0_dp * x**7.0_dp - 72072.0_dp * x**5.0_dp + 27720.0_dp * x**3.0_dp - 2520.0_dp * x) / 128.0_dp
    case(9)
        Val = (109395.0_dp * x**8.0_dp - 180180.0_dp * x**6.0_dp + 90090.0_dp * x**4.0_dp - 13860.0_dp * x**2.0_dp + 315.0_dp) / 128.0_dp
    case(10)
        Val = (461890.0_dp * x**9.0_dp - 875160.0_dp * x**7.0_dp + 540540.0_dp * x**5.0_dp - 120120.0_dp * x**3.0_dp + 6930.0_dp * x) / 256.0_dp
    end select


end function dLegendre_dx

real(dp) function SE_BF(Spectral, PolyOrder, N, Xi) result(SB)

    type(t_spectral)    :: Spectral
    integer             :: PolyOrder
    integer             :: N
    real(dp)            :: Xi

    SB = 0.0_dp
    if (abs(Xi - Spectral%Xi(N)) < 1.0E-12_dp) then
        SB = 1.0_dp
    else
        SB = 1.0_dp / (int(PolyOrder,dp) * (int(PolyOrder,dp) + 1) * Legendre(PolyOrder, Spectral%Xi(N))) * &
        Lobatto1D(PolyOrder-1, Xi) * (Xi ** 2.0_dp - 1.0_dp) / (Xi - Spectral%Xi(N)) 
    end if

end function SE_BF

real(dp) function SE_dBF(Spectral, PolyOrder, N, Xi) result(SB_D)

    type(t_spectral)    :: Spectral
    integer             :: PolyOrder
    integer             :: N, p, q
    real(dp)            :: Xi, TempValue

    SB_D = 0.0_dp
    do p = 1, 1 + PolyOrder
        if (p /= N) then
            TempValue = 1.0_dp / (Spectral%Xi(N) - Spectral%Xi(p))
            do q = 1, 1 + PolyOrder
                if (q /= N .and. q /= p) then
                    TempValue = TempValue * (xi - Spectral%Xi(q)) / (Spectral%Xi(N) - Spectral%Xi(q))
                end if
            end do
            SB_D = SB_D + TempValue
        end if
    end do

end function SE_dBF

subroutine Calculate_SE_Vec_Mat(&
    SE, Quad, QuadBound, SF_vec, dSF_dxi, dSF_deta, SF_vec_Bound, dSF_vec_Bound, SFMatrix, Poly)

    type(t_spectral)            :: SE
    type(t_Quadrature)          :: Quad, QuadBound
    real(dp), dimension(:,:)    :: SF_vec_Bound, dSF_vec_Bound
    real(dp), dimension(:,:)    :: SF_vec, dSF_dxi, dSF_deta
    real(dp), dimension(:,:,:)  :: SFMatrix
    integer                     :: Poly
    integer                     :: NBASIS, i, j, k, n

    NBASIS = Poly + 1

    SF_vec = 0.0_dp
    dSF_dxi = 0.0_dp
    dSF_deta = 0.0_dp
    SFMatrix = 0.0_dp
    SF_vec_Bound = 0.0_dp
    dSF_vec_Bound = 0.0_dp

    do i = 1, Quad%NoPoints
        n = 1
        do j = 1, NBASIS
            do k = 1, NBASIS
                SF_vec(i,n) = SE_BF(SE, Poly, j, Quad%Xi(i)) * SE_BF(SE, Poly, k, Quad%Eta(i))
                dSF_dxi(i,n) = SE_dBF(SE, Poly, j, Quad%Xi(i)) * SE_BF(SE, Poly, k, Quad%Eta(i))
                dSF_deta(i,n) = SE_BF(SE, Poly, j, Quad%Xi(i)) * SE_dBF(SE, Poly, k, Quad%Eta(i))
                n = n + 1
            end do
        end do

        SFMatrix(i,:,:) = outer_product(SF_vec(i,:), SF_vec(i,:))
    end do

    do i = 1, QuadBound%NoPoints
        n = 1
        do j = 1, NBASIS
            SF_vec_Bound(i,n) = SE_BF(SE,Poly,j,QuadBound%Xi(i))
            dSF_vec_Bound(i,n) = SE_dBF(SE,Poly,j,QuadBound%Xi(i))
            n = n + 1
        end do
    end do


end subroutine Calculate_SE_Vec_Mat

subroutine Calculate_SE_Vol_Basis(SE, Quad, SF_vec, SF_mat, dSF_dxi, dSF_deta, Poly)

    type(t_spectral)            :: SE
    type(t_Quadrature)          :: Quad
    real(dp), dimension(:,:)    :: SF_vec, dSF_dxi, dSF_deta
    real(dp), dimension(:,:,:)  :: SF_mat
    integer                     :: Poly

    integer                     :: NBASIS, i, j, k, n

    NBASIS = Poly + 1
    SF_vec = 0.0_dp
    dSF_dxi = 0.0_dp
    dSF_deta = 0.0_dp
    SF_mat = 0.0_dp

    do i = 1, Quad%NoPoints
        n = 1
        do j = 1, NBASIS
            do k = 1, NBASIS
                SF_vec(i,n) = SE_BF(SE, Poly, j, Quad%Xi(i)) * SE_BF(SE, Poly, k, Quad%Eta(i))
                dSF_dxi(i,n) = SE_dBF(SE, Poly, j, Quad%Xi(i)) * SE_BF(SE, Poly, k, Quad%Eta(i))
                dSF_deta(i,n) = SE_BF(SE, Poly, j, Quad%Xi(i)) * SE_dBF(SE, Poly, k, Quad%Eta(i))
                n = n + 1
            end do
        end do

        SF_mat(i,:,:) = outer_product(SF_vec(i,:), SF_vec(i,:))
    end do

end subroutine Calculate_SE_Vol_Basis

subroutine Calculate_SE_Vec_Mat_3D(&
    SE, Quad, Qsurf, SF_vec, dSF_dxi, dSF_deta, dSF_dzeta, SFMatrix, &
    SF_vec_bound, dSF_dxi_bound, dSF_deta_bound, Poly)

    type(t_spectral)                        :: SE
    type(t_Quadrature)                      :: Quad, Qsurf
    real(dp), dimension(:,:)                :: SF_vec, dSF_dxi, dSF_deta, dSF_dzeta
    real(dp), dimension(:,:), allocatable   :: SF_vec_bound, dSF_dxi_bound, dSF_deta_bound
    real(dp), dimension(:,:,:)              :: SFMatrix
    integer                                 :: Poly

    integer                                 :: NBASIS, gp, i, j, k, n

    NBASIS = Poly + 1
    
    SF_vec = 0.0_dp
    dSF_dxi = 0.0_dp
    dSF_deta = 0.0_dp
    SFMatrix = 0.0_dp

    SF_vec_bound = 0.0_dp
    dSF_dxi_bound = 0.0_dp
    dSF_deta_bound = 0.0_dp

    !$OMP PARALLEL DO PRIVATE(gp,n,i,j,k) SHARED(SF_vec,dSF_dxi,dSF_deta,SFMatrix,SE,Quad,Poly,NBASIS)
    do gp = 1, Quad%NoPoints
        n = 1
        do i = 1, NBASIS
            do j = 1, NBASIS
                do k = 1, NBASIS
                    SF_vec(gp,n) = SE_BF(SE,Poly,i,Quad%Xi(gp)) * SE_BF(SE,Poly,j,Quad%Eta(gp)) * SE_BF(SE,Poly,k,Quad%Zeta(gp))
                    dSF_dxi(gp,n) = SE_dBF(SE,Poly,i,Quad%Xi(gp)) * SE_BF(SE,Poly,j,Quad%Eta(gp)) * SE_BF(SE,Poly,k,Quad%Zeta(gp))
                    dSF_deta(gp,n) = SE_BF(SE,Poly,i,Quad%Xi(gp)) * SE_dBF(SE,Poly,j,Quad%Eta(gp)) * SE_BF(SE,Poly,k,Quad%Zeta(gp))
                    dSF_dzeta(gp,n) = SE_BF(SE,Poly,i,Quad%Xi(gp)) * SE_BF(SE,Poly,j,Quad%Eta(gp)) * SE_dBF(SE,Poly,k,Quad%Zeta(gp))
                    n = n + 1
                end do
            end do
        end do

        SFMatrix(gp,:,:) = outer_product(SF_vec(gp,:), SF_vec(gp,:))
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(gp,n,i,j) SHARED(SF_vec_bound,dSF_dxi_bound,dSF_deta_bound,SE,Qsurf,Poly,NBASIS)
    do gp = 1, Qsurf%NoPoints
        n = 1
        do i = 1, NBASIS
            do j = 1, NBASIS
                SF_vec_bound(gp,n) = SE_BF(SE, Poly, i, Qsurf%Xi(gp)) * SE_BF(SE, Poly, j, Qsurf%Eta(gp))
                dSF_dxi_bound(gp,n) = SE_dBF(SE, Poly, i, Qsurf%Xi(gp)) * SE_BF(SE, Poly, j, Qsurf%Eta(gp))
                dSF_deta_bound(gp,n) = SE_BF(SE, Poly, i, Qsurf%Xi(gp)) * SE_dBF(SE, Poly, j, Qsurf%Eta(gp))
                n = n + 1
            end do
        end do
    end do
    !$OMP END PARALLEL DO

end subroutine Calculate_SE_Vec_Mat_3D
end module m_spectral