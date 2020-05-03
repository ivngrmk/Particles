program MonteKarlo
    !То, что нужно для работы графики.
    use dflib
    type (wxycoord) wxy
    integer(2) logic
    !Для получения случаных чисел из файла.
    integer, parameter :: LastRandomNumber = 2400
    integer CurrentRandomNumber
    !Параметры системы.
    real, parameter :: av_t = 1.0 !Время свободного свободного пробега. ЗДЕСЬ СОВПАДАЕТ С ДЛИННОЙ, ТАК КАК СКОРОСТЬ ЧАСТИЦ = 1 У.Е. ИНАЧЕ РАБОТАТЬ БУДЕТ НЕКОРРЕКТНО.
    real, parameter :: pi = 3.14159265358979
    real t, t_sim
    real D !Коэффициент диффуции
    !Параметры симуляции
    integer, parameter :: NumberOfParticles = 100
    real, parameter :: MaxT = 100.0 !Время в у.е. в течение которого проводится симуляция.
    !Счётчики.
    integer i, j
    real(16) l_sum
    integer(8) count
    !Создание типа ЧАСТИЦА.
    type particle
        real x  !Координаты .
        real y
        real z
        real(16) vx !Направляющие косинусы.
        real(16) vy
        real(16) vz
    end type particle
    type(particle) part
    !Вспомогательные.
    real l, temp
    real tau
    real, parameter :: dt = 1.0
    !Нужно для получения случаных чисел.
    call RANDOM_SEED()
    !Начальные параметры системы
    D = av_t * 1.0 / 3.0    !!!Так как скорость частицы численно равно единице, то время пробега расстояния l численно равно этому расстоянию!!!
    !Создание графического окна и прорисовка в нём окна с графиком.
    call GraphicWindow()
    call GraphicAxes()
    !Осуществление получения случайных чисел из файла.
    CurrentRandomNumber = 1; open(2,file = "RandomNumbers.txt")
    !Открываем файл для записи необходимых данных.
    open(1, file= "out.txt")
    
    call MoveTo_w(DBLE(0.0), DBLE(0.0), wxy); logic = SetColor(4);

    call RANDOM_SEED()
    
    !Считаем дистанцию, которую i-ая частица прошло за время Tmax.
    t = 0
    !Начальное состояние этой частицы.
        
    !Считаем дистанцию, на которая одалится от начала координат i-ая частица за время MaxT.
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
    
    real function Get_Random_Time2Strike() !Получание случайного времени пробега частицы.
        real r
        r = GetRandomNumberFromFile()
        Get_Random_Time2Strike = -av_t * log( 1 - r )
    end function Get_Random_Time2Strike
    
    real function Get_Random_CosTetha()    !Получение угла в плоскости соударения частиц.
        real r
        r = GetRandomNumberFromFile()
        Get_Random_CosTetha = 1 - 2 * r
    end function Get_Random_CosTetha
    
    real function Get_Random_Phi()      !Получаение угла в плоскости перпенддикулярной линии соударения частиц.
        real r
        r = GetRandomNumberFromFile()
        Get_Random_Phi = 2 * pi * r
    end function Get_Random_Phi
    
    type(particle) function ConvertPart( old_part, cos_tetha,  phi ) !Построение нового вектора скорости по углам соударения.
        type(particle) old_part
        real cos_tetha, phi
        ConvertPart.x = old_part.x; ConvertPart.y = old_part.y; ConvertPart.z = old_part.z;
        ConvertPart.vx = cos_tetha * old_part.vx - ( old_part.vy * sin(phi) - old_part.vx * old_part.vz * cos(phi) ) * sqrt ( (1 - (cos_tetha)**2) / (1 - (old_part.vz)**2) )
        ConvertPart.vy = cos_tetha * old_part.vy + ( old_part.vx * sin(phi) + old_part.vy * old_part.vz * cos(phi) ) * sqrt ( (1 - (cos_tetha)**2) / (1 - (old_part.vz)**2) )
        ConvertPart.vz = cos_tetha * old_part.vz - ( 1 - (old_part.vz)**2 ) * cos(phi) * sqrt ( (1 - (cos_tetha)**2) / (1 - (old_part.vz)**2) )
    end function ConvertPart
    
    subroutine GraphicWindow() !Создаёт окно для вывода графики. фон - белый, граница - чёрная, само окно - белое.
        logical(2) bool2
        integer Pxl, Pxr, Pyu, Pyd
        Pxl = 100; Pyl = 50; Pxr = 900; Pyr = 650
        bool2 = SetBkColor(15); call ClearScreen(0) !Окраска всего экрана в белый.
        bool2 = SetColor(0); bool2 = Rectangle($GBORDER,Pxl-1, Pyl-1, Pxr+1, Pyr+1) !Создание границ окна.
        call SetViewPort(Pxl, Pyl, Pxr, Pyr) !Создание рабочей области.
        bool2 = SetBkColor(15); call ClearScreen(1) !Окраска рабочей области под график в белый. 1 - ссылка на эту рабочую область.
    end subroutine GraphicWindow
    
    subroutine GraphicAxes()
        real xl, yl, xr, yr, scale_width
        real x, y
        xl = -0.1; yl = -0.1; xr = MaxT; yr = 70.0; scale_width = 0.1 !Обязательно должны содержать начало координат.
        bool2 = SetWindow(.TRUE., DBLE(xl), DBLE(yl), DBLE(xr), DBLE(yr))
        x = xl
        do while (ceiling(x) <= floor(xr)) !Градуировка шкалы абсцисс.
            call MoveTo_w(DBLE(ceiling(x)), DBLE(0.0 - scale_width), wxy)
            bool2 = LineTo_w(DBLE(ceiling(x)), DBLE(0.0 + scale_width))
            x = x + 1
        end do
        y = yl
        do while (ceiling(y) <= floor(yr)) !Градуировка шкалы ординат.
            call MoveTo_w(DBLE(0.0 - scale_width), DBLE(ceiling(y)), wxy)
            bool2 = LineTo_w(DBLE(0.0 + scale_width), DBLE(ceiling(y)))
            y = y + 1
        end do
        bool2 = SetColor(4) !Рисуем сами оси.
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