program MonteKarlo
    !��, ��� ����� ��� ������ �������.
    use dflib
    type (wxycoord) wxy
    integer(2) logic
    !��� ��������� �������� ����� �� �����.
    integer, parameter :: LastRandomNumber = 2400
    integer CurrentRandomNumber
    !��������� �������.
    real, parameter :: av_t = 1.0 !����� ���������� ���������� �������. ����� ��������� � �������, ��� ��� �������� ������ = 1 �.�. ����� �������� ����� �����������.
    real, parameter :: pi = 3.14159265358979
    real t, t_sim
    real D !����������� ��������
    !��������� ���������
    integer, parameter :: NumberOfParticles = 100
    real, parameter :: MaxT = 100.0 !����� � �.�. � ������� �������� ���������� ���������.
    !��������.
    integer i, j
    real(16) l_sum
    integer(8) count
    !�������� ���� �������.
    type particle
        real x  !���������� .
        real y
        real z
        real(16) vx !������������ ��������.
        real(16) vy
        real(16) vz
    end type particle
    type(particle) part
    !���������������.
    real l, temp
    real tau
    real, parameter :: dt = 1.0
    !����� ��� ��������� �������� �����.
    call RANDOM_SEED()
    !��������� ��������� �������
    D = av_t * 1.0 / 3.0    !!!��� ��� �������� ������� �������� ����� �������, �� ����� ������� ���������� l �������� ����� ����� ����������!!!
    !�������� ������������ ���� � ���������� � �� ���� � ��������.
    call GraphicWindow()
    call GraphicAxes()
    !������������� ��������� ��������� ����� �� �����.
    CurrentRandomNumber = 1; open(2,file = "RandomNumbers.txt")
    !��������� ���� ��� ������ ����������� ������.
    open(1, file= "out.txt")
    
    call MoveTo_w(DBLE(0.0), DBLE(0.0), wxy); logic = SetColor(4);

    call RANDOM_SEED()
    
    !������� ���������, ������� i-�� ������� ������ �� ����� Tmax.
    t = 0
    !��������� ��������� ���� �������.
        
    !������� ���������, �� ������� �������� �� ������ ��������� i-�� ������� �� ����� MaxT.
    do while (t < MaxT)
        l_sum = 0
        do i = 1,NumberOfParticles
            t_sim = 0
            part.x = 0.0; part.y = 0.0; part.z = 0.0;
            part.vx = 1.0; part.vy = 0.0; part.vz = 0.0;
            count = 0
            do while (t_sim < t)
                count = count + 1
                tau = Get_Random_Time2Strike()
                if (t_sim + tau < t) then
                    part.x = part.x + part.vx * tau
                    part.y = part.y + part.vy * tau
                    part.z = part.z + part.vz * tau
                    part = ConvertPart( part, Get_Random_CosTetha(), Get_Random_Phi() )
                else
                    tau = (t - t_sim)
                    part.x = part.x + part.vx * tau
                    part.y = part.y + part.vy * tau
                    part.z = part.z + part.vz * tau
                end if
                t_sim = t_sim + tau
            end do
            l = sqrt(part.x**2 + part.y**2 + part.z**2)
            l_sum = l_sum + l
        end do
        l = l_sum / NumberOfParticles
        logic = LineTo_w( DBLE(t), DBLE(l) )
        write(1,*) t, " ", l
        t = t + dt
    end do
    
    t = 0
    call MoveTo_w(DBLE(0.0), DBLE(0.0), wxy); logic = SetColor(1);
    do while (t < MaxT)
        t = t + dt
        logic = LineTo_w( DBLE(t), DBLE( sqrt(6*D*t) ) )
    end do

    close(1)
    close(2)
    
    contains
    
    real function Get_Random_Time2Strike() !��������� ���������� ������� ������� �������.
        real r
        r = GetRandomNumberFromFile()
        Get_Random_Time2Strike = -av_t * log( 1 - r )
    end function Get_Random_Time2Strike
    
    real function Get_Random_CosTetha()    !��������� ���� � ��������� ���������� ������.
        real r
        r = GetRandomNumberFromFile()
        Get_Random_CosTetha = 1 - 2 * r
    end function Get_Random_CosTetha
    
    real function Get_Random_Phi()      !���������� ���� � ��������� ����������������� ����� ���������� ������.
        real r
        r = GetRandomNumberFromFile()
        Get_Random_Phi = 2 * pi * r
    end function Get_Random_Phi
    
    type(particle) function ConvertPart( old_part, cos_tetha,  phi ) !���������� ������ ������� �������� �� ����� ����������.
        type(particle) old_part
        real cos_tetha, phi
        ConvertPart.x = old_part.x; ConvertPart.y = old_part.y; ConvertPart.z = old_part.z;
        ConvertPart.vx = cos_tetha * old_part.vx - ( old_part.vy * sin(phi) - old_part.vx * old_part.vz * cos(phi) ) * sqrt ( (1 - (cos_tetha)**2) / (1 - (old_part.vz)**2) )
        ConvertPart.vy = cos_tetha * old_part.vy + ( old_part.vx * sin(phi) + old_part.vy * old_part.vz * cos(phi) ) * sqrt ( (1 - (cos_tetha)**2) / (1 - (old_part.vz)**2) )
        ConvertPart.vz = cos_tetha * old_part.vz - ( 1 - (old_part.vz)**2 ) * cos(phi) * sqrt ( (1 - (cos_tetha)**2) / (1 - (old_part.vz)**2) )
    end function ConvertPart
    
    subroutine GraphicWindow() !������ ���� ��� ������ �������. ��� - �����, ������� - ������, ���� ���� - �����.
        logical(2) bool2
        integer Pxl, Pxr, Pyu, Pyd
        Pxl = 100; Pyl = 50; Pxr = 900; Pyr = 650
        bool2 = SetBkColor(15); call ClearScreen(0) !������� ����� ������ � �����.
        bool2 = SetColor(0); bool2 = Rectangle($GBORDER,Pxl-1, Pyl-1, Pxr+1, Pyr+1) !�������� ������ ����.
        call SetViewPort(Pxl, Pyl, Pxr, Pyr) !�������� ������� �������.
        bool2 = SetBkColor(15); call ClearScreen(1) !������� ������� ������� ��� ������ � �����. 1 - ������ �� ��� ������� �������.
    end subroutine GraphicWindow
    
    subroutine GraphicAxes()
        real xl, yl, xr, yr, scale_width
        real x, y
        xl = -0.1; yl = -0.1; xr = MaxT; yr = 70.0; scale_width = 0.1 !����������� ������ ��������� ������ ���������.
        bool2 = SetWindow(.TRUE., DBLE(xl), DBLE(yl), DBLE(xr), DBLE(yr))
        x = xl
        do while (ceiling(x) <= floor(xr)) !����������� ����� �������.
            call MoveTo_w(DBLE(ceiling(x)), DBLE(0.0 - scale_width), wxy)
            bool2 = LineTo_w(DBLE(ceiling(x)), DBLE(0.0 + scale_width))
            x = x + 1
        end do
        y = yl
        do while (ceiling(y) <= floor(yr)) !����������� ����� �������.
            call MoveTo_w(DBLE(0.0 - scale_width), DBLE(ceiling(y)), wxy)
            bool2 = LineTo_w(DBLE(0.0 + scale_width), DBLE(ceiling(y)))
            y = y + 1
        end do
        bool2 = SetColor(4) !������ ���� ���.
        call MoveTo_w(DBLE(xl), DBLE(0.0), wxy)
        bool2 = LineTo_w(DBLE(xr), DBLE(0.0))
        call MoveTo_w(DBLE(0.0), DBLE(yl), wxy)
        bool2 = LineTo_w(DBLE(0.0), DBLE(yr))
    end subroutine GraphicAxes
    
    real function GetRandomNumberFromFile()
        real r
        read(2,*) r
        CurrentRandomNumber = CurrentRandomNumber + 1
        if (CurrentRandomNumber > LastRandomNumber) then
            CurrentRandomNumber = 1
            close(2)
            open(2,file = "RandomNumbers.txt")
        end if
        GetRandomNumberFromFile = r
    end function GetRandomNumberFromFile
    
end