\documentclass[12pt]{report}
\usepackage[a4paper, left=12mm, right=10mm, top=13mm, bottom=15mm]{geometry}
\usepackage[unicode]{hyperref}
\usepackage[utf8]{inputenc}
\usepackage[russian]{babel}
\usepackage[english]{babel}
\usepackage{indentfirst}
\usepackage{minted}
\usepackage{xcolor}
\usepackage{keyval}
\usepackage{fancyvrb}
\usepackage{float}
\usepackage{ifthen}
\usepackage{calc}

\documentclass[a4paper, 14pt]{extreport}
\usepackage{scrextend}
\usepackage{amsmath}
\usepackage[a4paper, left = 20mm, right = 20mm, top = 13mm, bottom = 15mm]{geometry}
\usepackage[unicode]{hyperref}
\usepackage{indentfirst}
\usepackage{minted}
\usepackage{xcolor}
\usepackage{keyval}
\usepackage{fancyvrb}
\usepackage{float}
\usepackage{ifthen}
\usepackage{calc}
\usepackage{ifplatform}
\AtBeginEnvironment{minted}{\singlespacing%
    \fontsize{14}{5}\selectfont}
\usepackage{listings}
\usepackage{fancyvrb}
\usepackage{setspace,amsmath}
\urlstyle{same}
\usepackage{graphicx}
\graphicspath{./}
\usepackage[utf8]{inputenc}
\usepackage[russian]{babel}
\usepackage{pgfplots}
\pgfplotsset{compat=1.9}
\usepackage[usenames]{color}
\usepackage{colortbl}
\usepackage{ifplatform}
\usepackage{etoolbox}
\AtBeginEnvironment{minted}{\singlespacing%
    \fontsize{13}{5}\selectfont}
\usepackage{listings}
\usepackage{fancyvrb}
\usepackage{setspace,amsmath}


\title{{\textbf{Введение в Численные Методы\\ Аналитический отчёт по практическому заданию}}}
\author{Выполнила студентка 208 группы ВМК МГУ\\ Мазур Анастасия Вадимовна}
\date{}

\begin{document}
\maketitle

\begin{center}
    Var 2 ex 1
\end{center}

\section*{Математическая постановка задачи}

Функция $f(x)$  задана таблично на отрезке [0, a] в точках $x_{i}$, $x_{i} = ih$, \\
$i = 0, 1, ..., n, h = a / n$

\begin{itemize}
\item{Построить интерполяционный многочлен по точкам $x$.} 
\item{Приблизить функцию по методу наименьших квадратов  полиномом   заданной степени $n$, $n < 9$. Оценить погрешность.} 
\item{Результаты сравнить.}
\end{itemize}

\noindent{Отрезок $[0, 2], \;\; n = 4$}\\
Таблица значений функции в точках: \\ \\
\indent\begin{tabular}{|c|c|c|} 
    \hline\( i \) & \( x \) & \( f(x) \) \\ 
    \hline 0 & 0 & 0 \\ 
    \hline 1 & \( 0.2 \) & \( 0.006732 \) \\ 
    \hline 2 & \( 0.4 \) & \( 0.058195 \) \\ 
    \hline 3 & \( 0.6 \) & \( 0.030482 \) \\ 
    \hline 4 & \( 0.8 \) & \( 0.387483 \) \\ 
    \hline 5 & \( 1 \) & \( 0.958924 \) \\ 
    \hline 6 & \( 1.2 \) & \( 0.48283 \) \\
    \hline 7 & \( 1.4 \) & \( 1.802771 \) \\ 
    \hline 8 & \( 1.6 \) & \( 4.052411 \) \\
    \hline 9 & \( 1.8 \) & \( 2.403475 \) \\
    \hline 10 & \( 2 \) & \( 4.352169 \) \\
    \hline
\end{tabular}

\section*{Используемые алгоритмы и формулы}
\subsection*{Построение интерполяционного многочлена \\ в форме Лагранжа}
\;\;\;\;\;\;Чтобы интерполировать функцию построим полином в форме Лагранжа. \\ Искомый полином $P_{n}(x)$ будет иметь следующий вид:

$$P_{n}(x) = \sum_{i = 0}^{n} f\left(x_{i}\right) Q_{n, i}(x),$$
где $Q_{n, i}(x)$ - полиномы степени $n$, "ориентированные" на точки $x_{i}$ \\ в том смысле, что
$$
Q_{n, i}(x)=\left\{\begin{array}{l}0, x=x_{j}, \; \; \forall j \neq i \\ 1, x=x_{i}\end{array}\right. $$

{Полиномы имееют вид: $$Q_{n, i}(x)=\prod_{j=0 \atop j \neq i}^{j=n} \frac{\left(x-x_{j}\right)}{\left(x_{i}-x_{j}\right)}$$}
или в нашем случае, когда $x_{i} = ih$, то есть известны значения в точках, расстояние между которыми фиксировано, то выражение можно упростить до следующей записи:
$$	
Q_{n, i}(x)=h^{-n} \cdot \prod_{j=0 \atop j \neq i}^{n} \frac{(x-j h)}{(i-j)}
$$

Учитывая, что полином в форме Лагранжа $P_{n}(x)$ представляет собой линейную комбинацию алгебраических уравнений $f(x_{i})Q_{n, i}(x), Q_{n, i}(x)$ - полиномы степени n, можно утверждать, что $P_{n}(x)$ будет иметь степень не более $n$.
\\ \\
\indent{Данные формулы будут далее использоваться в программной реализации.}
\subsection*{Приближение функции методом наименьших квадратов}

\;\;\;\;\;\;В методе наименьших квадратов аппроксимирующая функция $y(x)$ ищется в виде следующей суммы:
$$F(x)=\sum_{k=0}^{m} a_{k} \varphi_{k}(x), \; \; m<n$$

{В каждой точке сетки $x_{i}$ можно подсчитать погрешность:}
$$\delta_{i}=y_{i}-F\left(x_{i}\right)=y_{i}-\sum_{k=0}^{m} a_{k} \varphi_{k}\left(x_{i}\right), \; \; \; i=0,1,2, \ldots, n$$

Сумма квадратов этих величин называется суммарной квадратичной погрешностью
$$J=\sum_{i=0}^{n} \delta_{i}^{n}=\sum_{i=0}^{n}\left(y_{i}-\sum_{k=0}^{m} a_{k} \varphi_{k}\left(x_{i}\right)\right)^{2}
$$

Главной задачей является подобрать такие коэффициенты $a_{k}$, чтобы суммарная квадратичная погрешность была минимальной.
\\ \indent{Таким образом, построение наилучшего приближения сводится к классической задаче математического анализа об экстремуме функции нескольких переменных. Необходимым условием экстремума является равенство нулю в экстремальной точке всех первых частных производных функции.}

$$\frac{\partial J}{\partial a_{e}}=-2 \sum_{i=0}^{n}\left(y_{i} \sum_{k=0}^{m} a_{k} \varphi_{k}\left(x_{i}\right)\right) \varphi_{L}\left(y_{i}\right)=0, \; \; \; l = 0, 1, \ldots, m.$$

Оставим члены, содержащие $a_{k}$, слева и поменяем в них порядок суммирования по индексам $i$ и $k$. Члены, содержащие $y_{i}$, перенесем направо. В результате уравнения примут вид:
$$\sum_{k=0}^{m} \gamma_{l k} a_{k}=b_{l}, \; \; \; l=0,1, \ldots, m,$$ где
$$\begin{aligned} \gamma_{l k} &=\sum_{i=0}^{n} \varphi_{l}\left(x_{i}\right) \varphi_{k}\left(x_{i}\right) \\ b_{l} &=\sum_{i=0}^{n} \varphi_{l}(x_{i}) y_{i} \end{aligned}$$

Мы получили систему линейных алгебраических уравнений, в которой роль неизвестных играют искомые коэффициенты разложения $a_{0}, a_{1}, \ldots, a_{m}$. Используя найденные коэффициенты разложения, мы сможем построить наилучшее приближении сеточной функции по методу наименьших квадратов.
\\ \\ \indent{Данные формулы будут далее использоваться в программной реализации.}

\section*{Цифровое представление результатов}
\;\;\;\;\;\;Для того, чтобы произвести необходимые вычисления и проанализировать результаты, мною была написана программа. \\ 
\indentВ данном отчёте представлены лишь основные функции и результат работы программы на входных данных из условия, однако при желании с полным кодом можно ознакомиться отдельно и произвести тестирование на других входных данных (см. архив).\\
\\
\indent\textbf{Принцип работы программы} \\
\indent Известные точки записываются записываются в файл $dots.txt$, каждая точка на отдельной строке, абсцисса и ордината отделены пробельным символом. Точки, значения которых хотим узнать, записываются в файл $input.txt$. \\ \indent После запуска программы в файлах $output\_Lagrange.txt$ и \\ $output\_LeastSquares.txt$ будут записаны прогнозируемые значения, посчитанные по этим методам соотвественно. \\
\indentТакой подход позволяет изменять входные данные, не изменяя программного кода.
\\ \indentЗатем файлы можно сравнить с помощью следующей команды терминала: \textit{meld}, она наиболее наглядно показывает отличия в файлах.

\noindent\includegraphics[width =\textwidth]{meld}
Сравнение предсказанных значений в произвольных точках отрезка на основе вычислений по разным методам c помощью \textit{meld}\\
Слева - интерполяция полиномом Лагранжа \\
Справа - аппроксимация методом наименьших квадратов \\

\indent \textbf{Функции}

\vspace{-\medskipamount}
\vspace{-2\baselineskip}
\vspace{10pt}
\usemintedstyle[c]{friendly} 
\begin{minted}{c}
long double
y_Lagrange(long double *x, long double *y, int num, long double x0) {
    long double L;
    int i, j;
    long double P;

    L = 0;
    for (i = 0; i < num; i++) {
        P = 1;
        for (j = 0; j < num; j++) {
            if (i != j) {
                P *= (x0 - x[j]) / (x[i] - x[j]);              
            }
        }
        L += y[i] * P;
    }
    return L;
}
\end{minted}
\\
\vspace{-\medskipamount}
\vspace{-2\baselineskip}
\vspace{10pt}
\usemintedstyle[c]{friendly} 
\begin{minted}{c}
long double
y_LeastSquares(long double *x, long double* y, int num, long double x0) {
    int i, j, k;

    for (i = 0; i < n + 1; i++){
        a[i]=0;
        b[i]=0;
            for ( j = 0; j < n + 1; j++){
                sum[i][j] = 0;
            }
    }

    for (i = 0; i < n + 1; i++) {
        for (j = 0; j < n + 1; j++) {
            sum[i][j] = 0;
            for (k = 0; k < num; k++) {
                sum[i][j] += pow(x[k], i + j);
            }
        }
    }
    
    for (i = 0; i < n + 1; i++) {
        for (k = 0; k < num; k++) {
            b[i] += pow(x[k], i) * y[k];
        }
    }

    long double temp = 0;
    for (i = 0; i < n + 1; i++) {
        if (sum[i][i] == 0) {
            for (j = 0; j < n + 1; j++) {
                if (j == i) {
                    continue;
                }
                if (sum[j][i] != 0 && sum[i][j] != 0) {
                    for (k = 0; k < n + 1; k++) {
                        temp = sum[j][k];
                        sum[j][k] = sum[i][k];
                        sum[i][k] = temp;
                    }
                    temp = b[j];
                    b[j] = b[i];
                    b[i] = temp;
                    break;
                }
            }
        }
    }

    for (k = 0; k < n + 1; k++) {
        for (i = k + 1; i < n + 1; i++) {
            if (sum[k][k] == 0) {
                printf("\n Solution doesn't exist! \n");
                return 0;
            }
            long double M = sum[i][k] / sum[k][k];
            for (j = k; j < n + 1; j++) {
                sum[i][j] -= M * sum[k][j];
        }
        b[i] -= M * b[k];
    }
}

    for (i = n; i >= 0; i--) {
        long double s = 0;
        for (j = i; j < n + 1; j++){
            s += sum[i][j] * a[j];
        }
        a[i] = (b[i] - s) / sum[i][i];
    }

    long double result = 0;

    for (i = 0; i < n + 1; i++){
        result += a[i] * pow(x0, i);
    }
    return result;
}

\end{minted}

\section*{Графическое представление результатов}

\includegraphics[width =\textwidth]{graphic}
Синие точки - точки, известные из условия \\
Оранжевая кривая - кривая интерполирующего многочлена в \\ форме Лагранжа \\
Зелёная кривая - кривая, построенная по методу наименьших квадратов

\section*{Анализ результатов}
\;\;\;\;\;\;Ключевым отличием двух рассматриваемых методов является то, что полином Лагранжа интерполирует функцию $f(x)$, точки которой нам известны, а метод наименьших квадратов эту функцию аппроксимирует. \\
\indentТо есть для полинома Лагранжа важно, чтобы полученная интерполяционная функция строго проходила через известные узлы, однако вне известных точек функция сильно "скачет". Такой разброс значений приводит к тому, что полученная интерполяционная кривая плохо характеризует поведение исходной функции $f(x)$ в целом. Это усложняет прогнозирование значений и дальнейшую работу с функцией, затрудняет визуальное восприятие графика.\\
\indentМетод наименьших квадратов, напротив, на выходе даёт нам такую функцию, которая лишь приближает $f(x)$, то есть может совпадать с рассматриваемой функцией на очень маленьком наборе точек, но зато полученная кривая довольно удобна и наглядна. 

\section*{Дополнительно: применимость методов}
\;\;\;\;\;\;В данный момент метод наименьших квадратов активно применяется в анализе данных. Этот метод помогает аппроксимировать различные экспериментальные данные в разных областях, например, в экономике. Это очень сильный инструмент, потому что позволяет сделать довольное точное прогнозирование новых значений. \\
\indentЧто касается интерполяции с помощью полинома Лагранжа, то такой метод тоже пользуется популярностью в современном мире. Например, метод был использован для вычислений при разработке двигателей. Построенный полином помог сократить количество реальных тестирований.

\section*{Источники и ресурсы}
Вводные лекции по численным методам (Д.П. Костомаров, A.П. Фаворский) \\
Для построения графика использовался ресурc \textit{www.geogebra.com}
\end{document}


\title{ \textbf{Отчётная работа\\ по практикуму на ЭВМ}}
\author{Выполнила студентка 208 группы ВМК МГУ\\ Мазур Анастасия Вадимовна}
\date{}

\begin{document}
\maketitle

\definecolor{code}{HTML}{F4FAFC}
\definecolor{terminal}{HTML}{A9C3C8 }

\makeatletter
\def\verbatim{\footnotesize\@verbatim \frenchspacing\@vobeyspaces \@xverbatim}
\makeatother

\clearpage
\renewcommand{\contentsname}{ОГЛАВЛЕНИЕ}
\tableofcontents
\newpage

\chapter*{Задание 6}
\addcontentsline{toc}{chapter}{Задание 6}
\large
\textit{Проконтролировать, допускается ли инициализация переменных \\ при описании; происходит ли инициализация по умолчанию.}

\section*{Введение}
\addcontentsline{toc}{section}{Введение}

\text{\textbf{Инициализация переменной} - это процесс присваивания этой переменной} \newline \text{некоторого начального значения, которое может быть изменено в ходе программы.}
\newline \\
\indent \text{Несмотря на то, что во многих языках программирования описанным,} \newline 
\text{но явно не инициалазированным переменным по умолчанию присваивается} \newline \text{нулевое начальное значение, среди программистов всё же считается хорошим} \newline \text{тоном в явном виде инициализировать переменную. } \newline \indent \text{Такой стиль написания кода помогает избежать ошибок по невнимательности} \newline \text{или при работе в команде.}

\newline \newline

\section*{Постановка задачи}
\addcontentsline{toc}{section}{Постановка задачи}

\text{В процессе выполнения данной задачи необходимо опытным путём выяснить,} \newline \text{как работает инициализация переменных в языке Си:} \newline \text{- Можно ли инициализировать переменные сразу при их описании?} \newline \text{- Можно ли вовсе не инициализировать переменные? Если да, то определить,} \newline \text{какое начальное значение компилятор присвоит переменным по умолчанию.}
\newline \newline 
\indent \text{Примечание: для проверки программного кода, приведенного в листингах} \newline \text{использовалась следующая версия компилятора - \textbf{gcc version 9.3.0}} \newline

\newpage

\section*{Исследование задачи}
\addcontentsline{toc}{section}{Исследование задачи}
\text{Для того, чтобы ответить на вопросы задачи, отдельно рассмотрим переменные} \newline \text{различных классов и типов. Так мы сможем сделать необходимые выводы об}
\newline \text{особенностях их инициализации.}
\newline \newline
\indent\text{Скалярные переменные можно инициализировать в их определениях, помещая} \newline \text{после имени знак =}
\vspace{-\medskipamount}
\vspace{-2\baselineskip}
\vspace{25pt}
\usemintedstyle[c]{friendly} 
\begin{minted}[bgcolor=code]{c}
#include <stdio.h>

int 
main(void) {
    int a = 208;
    int a1 = 0xAB1; // число в шестнадцатиричной системе счисления
    float b = 2.5;
    double b1 = 4E-2; // экспоненциальная запись числа
    char c = '#';
    long g = 20*4; // вычисление значения, а затем преобразование типов
    long f = g+5L; // выражение использует ранее определённое значение
    
    printf("integer a = %i\n", a);
    printf("integer a1 = %i\n", a1);
    printf("float b = %f\n", b);
    printf("double b1 = %lf\n", b1);
    printf("character c = %c\n", c);
    printf("long integer g = %li\n", g);
    printf("long integer f = %li\n", f);
    
    return 0;
}
\end{minted}
\indent \text{Вывод:}

\vspace{-\medskipamount}
\vspace{-2\baselineskip}
\vspace{25pt}

\usemintedstyle[bash]{default} 
\begin{minted}[bgcolor=terminal]{bash}
integer a = 208
integer a1 = 2737
float b = 2.500000
double b1 = 0.040000
character c = #
long integer e = 80
long integer f = 85
\end{minted}
\vspace{-\medskipamount}
\vspace{-2\baselineskip}
\vspace{25pt}

\newline
\indent\text{Для внешних и статических переменных инициализирующие выражения} \newline \text{должны быть константными, при этом инициализация осуществляется только} \newline \text{один раз до начала выполнения программы.}
\newline 
\indent \text{Например, если приведённую выше программу дополнить следующим описанием} \newline \text{переменной: } 
\vspace{-\medskipamount}
\vspace{-2\baselineskip}
\vspace{25pt}
\usemintedstyle[c]{friendly} 
\begin{minted}[bgcolor=code]{c}
static int m = a1 - a;
\end{minted}

\text{то компилятор выдаст ошибку:}
\vspace{-\medskipamount}
\vspace{-2\baselineskip}
\vspace{25pt}
\usemintedstyle[bash]{default} 
\begin{minted}[bgcolor=terminal]{bash}
error: initializer element is not a compile-time constant
\end{minted}

\newpage
\text{Инициализация автоматических и регистровых переменных выполняется каждый} \newline \text{раз при входе в функцию или блок. Для таких переменных инициализирующее} \newline \text{выражение - не обязательно константное. Это может быть любое выражение,} \newline \text{использующее ранее определенные значения, включая даже и вызовы функций.}
\newline \indent \text{Рассмотрим пример:}
\vspace{-\medskipamount}
\vspace{-2\baselineskip}
\vspace{25pt}
\usemintedstyle[c]{friendly} 
\begin{minted}[bgcolor=code]{c}
#include <stdio.h>

void
func(void) {
    int a = 1;
    printf("a = %d\n", a++); // а == 2
    // но при следующем вызове функции a снова будет инициализирована 1;
}

int
main(void) {
    int i;
    for (i = 0; i < 3; i++) {
        func();
    }
    return 0;
}
\end{minted}
\indent \text{Вывод:}
\inputminted{cpp}{\jobname.cpp}

\vspace{-\medskipamount}
\vspace{-2\baselineskip}
\vspace{25pt}

\inputminted{cpp}{\jobname.cpp}

\usemintedstyle[bash]{default} 
\begin{minted}[bgcolor=terminal]{bash}
a = 1
a = 1
a = 1
\end{minted}
\indent \text{Что касается массивов, то их тоже можно инициализировать при описании.} \newline \text{Массив инициализируется в своём определении с помощью заключённого в фигурные} \newline \text{скобки списка инициализаторов, разделённых запятыми.}
\vspace{-\medskipamount}
\vspace{-2\baselineskip}
\vspace{25pt}
\usemintedstyle[c]{friendly} 
\begin{minted}[bgcolor=code]{c}
int arr[5] = {1,2,3,4,5};
\end{minted}
\indent \text{Если размер массива не указан, то длину массива компилятор вычисляет по числу} \newline \text{заданных инициализаторов.} \newline \indent\text{Символьные массивы можно инициализировать не используя конструкцию с} \newline \text{фигурными скобками, а задать значение с помощью строки.} \newline \indent \text{Эти записи эквивалентны:}
\vspace{-\medskipamount}
\vspace{-2\baselineskip}
\vspace{25pt}
\usemintedstyle[c]{friendly} 
\begin{minted}[bgcolor=code]{c}
char arr[] = "praktikum";
\end{minted}
\vspace{-\medskipamount}
\vspace{-2\baselineskip}
\vspace{8pt}
\usemintedstyle[c]{friendly} 
\begin{minted}[bgcolor=code]{c}
char arr[] = {'p','r','a','k','t','i','k','u','m','\0'};
// Си строки оканчиваются символом '\0'
\end{minted}
\newline \indent \text{Таким образом, мы рассмотрели примеры инициализации разнообразных типов} \newline \text{переменных и с уверенностью можем сделать вывод о том, что \underline{в языке Си можно}} \newline \text{\underline{инициализировать переменные сразу при их описании}, это упрощает код. Однако} \newline \text{необходимо учитывать особенности типа и производить только корректные} \newline \text{присваивания.}
\newpage
\indent\text{Что же происходит в тех случаях, когда переменная не инициализирована?} \newline \text{Произойдёт ли инициализация по умолчанию? Какое значение будет присвоено?}
\newline\text{Для проверки модифицируем написанную ранее программу: }
\vspace{-\medskipamount}
\vspace{-2\baselineskip}
\vspace{25pt}
\usemintedstyle[c]{friendly} 
\begin{minted}[bgcolor=code]{c}
#include <stdio.h>

int
main(void) {
    int a;      
    float b;      
    double b1;                                        
    char c;      
    long g;                                                            

    printf("integer a = %i\n", a);
    printf("float b = %f\n", b);
    printf("double b1 = %lf\n", b1);
    printf("character c = %c\n", c);
    printf("long integer g = %li\n", g);

    return 0;
}
\end{minted}
\indent \text{Компилируя с ключом -Wall получаем большое количество предупреждений:}
\inputminted{cpp}{\jobname.cpp}
\vspace{-\medskipamount}
\vspace{-2\baselineskip}
\vspace{12pt}

\inputminted{cpp}{\jobname.cpp}
\AtBeginEnvironment{minted}{\singlespacing%
    \fontsize{12}{5}\selectfont}
\usemintedstyle[bash]{default} 
\begin{minted}[bgcolor=terminal]{bash}
prog.c:11:32: warning: variable 'a' is uninitialized when used here [-Wuninitialized]
    printf("integer a = %i\n", a);
                               ^
prog.c:5:10: note: initialize the variable 'a' to silence this warning
    int a;
         ^
          = 0
prog.c:12:30: warning: variable 'b' is uninitialized when used here [-Wuninitialized]
    printf("float b = %f\n", b);
                             ^
prog.c:6:12: note: initialize the variable 'b' to silence this warning
    float b;
           ^
            = 0.0
prog.c:13:33: warning: variable 'b1' is uninitialized when used here [-Wuninitialized]
    printf("double b1 = %lf\n", b1);
                                ^~
prog.c:7:14: note: initialize the variable 'b1' to silence this warning
    double b1;
             ^
              = 0.0
prog.c:14:34: warning: variable 'c' is uninitialized when used here [-Wuninitialized]
    printf("character c = %c\n", c);
                                 ^
prog.c:8:11: note: initialize the variable 'c' to silence this warning
    char c;
          ^
           = '\0'
prog.c:15:38: warning: variable 'g' is uninitialized when used here [-Wuninitialized]
    printf("long integer g = %li\n", g);
                                     ^
prog.c:9:11: note: initialize the variable 'g' to silence this warning
    long g;
          ^
           = 0
\end{minted}
\AtBeginEnvironment{minted}{\singlespacing%
    \fontsize{14}{5}\selectfont}
\indent\text{Компилятор предлагает нам проинициализировать переменные нулём, но мы} \newline \text{не слушаем его и запускаем программу. Получаем следующий вывод в консоль: }
\inputminted{cpp}{\jobname.cpp}

\vspace{-\medskipamount}
\vspace{-2\baselineskip}
\vspace{2pt}

\inputminted{cpp}{\jobname.cpp}

\usemintedstyle[bash]{default} 
\begin{minted}[bgcolor=terminal]{bash}
integer a = 325111845
float b = 0.000000
double b1 = 0.000000
character c = 
long integer g = 0
\end{minted}
 \indent \text{Но не спешим радоваться полученному результату. Запустим программу ещё раз:}
\inputminted{cpp}{\jobname.cpp}

\vspace{-\medskipamount}
\vspace{-2\baselineskip}
\vspace{2pt}

\inputminted{cpp}{\jobname.cpp}

\usemintedstyle[bash]{default} 
\begin{minted}[bgcolor=terminal]{bash}
integer a = 226791461
float b = 0.000000
double b1 = 0.000000
character c = 
long integer g = 0
\end{minted}
\indent\text{Видим, что значение переменной типа integer изменилось. Оно будет разным при} \newline \text{каждом запуске программы. Таким образом, мы готовы ответить на поставленный} \newline \text{вопрос.}
\newline \indent \text {На этапе трансляции программы компилятор замечает, что переменные} \newline \text{не инициализированы и предупреждает программиста об этом. Однако это не} \newline \text{является фатальной ошибкой и трансляция не прерывается. Переменным типов} \newline \text{float, double, char, long и т.п. (в том числе и массивы этих типов)} \newline \text{присваивается нулевое значение. В переменной типа int же хранится, так скажем,} \newline \text{мусор - случайное значение, которое осталось в той области памяти, которая была} \newline \text{выделена под переменную. Это связано, в первую очередь, с оптимизацией работы} \newline \text{программы: если нет необходимости в инициализации, то незачем тратить ресурсы} \newline \text{для записи нулей.}
\newline\indent\text{Следовательно, делаем следующий вывод: во многих случаях компилятор поможет} \newline \text{программисту и присвоит нули (хотя это может не всегда быть уместно, например,} \newline \text{в теле программы переменная используется для умножения => результат будет}\newline\text{занулён, это бесполезно).} \newline \indent \text{\underline{Поэтому необходимо всегда инициализировать переменные исходя из поставленной}}\newline\text{\underline{задачи, чтобы не допускать ошибок.}}

\newpage
\section*{Литература и веб-ресурсы}
\addcontentsline{toc}{section}{Литература и веб-ресурсы}
\begin{itemize}
\item[$-$]
\text{Б. Керниган, Д. Ритчи} \newline \textbf{Язык программирования Си}
\item[$-$]
\text{Н.В.Вдовикина, И.В.Машечкин, А.Н.Терехин, В.В.Тюляева} \newline \textbf{Программирование в ОС UNIX на языке Си}
\item[$-$] \textbf{Онлайн справочник программиста на С и C++}
\end{itemize}
\newpage
\chapter*{Задание 13}
\addcontentsline{toc}{chapter}{Задание 13}

\large
\textit{Определите, каким образом расширяются unsigned- варианты "малых" целых \\ (до unsigned int как в "классическом" Си или сначала делается попытка \\ расширить до int, и только в случае неудачи - до} \textit{unsigned int).}

\section*{Введение}
\addcontentsline{toc}{section}{Введение}
\text{Сложно не согласиться с утверждением, что для того, чтобы быть хорошим} \newline \text{программистом, нужно понимать, каким образом выполняются процессы в}\newline \text{компьютере: какие арифметические операции и преобразования стоят между}\newline \text{написанным кодом и полученным результатом.}
\newline \indent\text{Именно поэтому в процессе изучения языка программирования Си особенно}\newline\text{важно уделить должное внимание преобразованию типов, которое, как оказалось,}\newline\text{может зависеть от реализации.}

\section*{Постановка задачи}
\addcontentsline{toc}{section}{Постановка задачи}
\indent\text{В процессе выполнения данного задания необходимо провести сравнительную}\newline\text{характеристику реализаций компиляторов языка Си в прошлом, выяснить, каким}\newline\text{образом работют компиляторы сейчас. А конкретнее, каким образом происходит}\newline\text{расширение малых беззнаковых целых чисел.}
\newpage
\section*{Исследование задачи}
\addcontentsline{toc}{section}{Исследование задачи}
\newline\text{Заглянем в историю: раньше не было единой спецификации, и поэтому разные}\newline \text{компиляторы поступали разными способами, то есть они отличались в своей}\newline \text{реализации.
Существует два основных типа реализаций: } \newline \text{\textbf{unsigned preserving (сохранение без знака)} и \textbf{value preserving (сохранение}} \newline\text{\textbf{значения).}}
\newline\indent\text{ANSI C использует для арифметических преобразований правила обработки с} \newline \text{сохранением значения, в то время как реализации K&R C склоняются к использованию} \newline \text{правил операций без знака. Например, следующий код по стандарту ANSI C делает} \newline \text{знаковое деление, где реализация unsigned-preserving делала бы деление без знака:}
\vspace{-\medskipamount}
\vspace{-2\baselineskip}
\vspace{2pt}
\usemintedstyle[c]{friendly} 
\begin{minted}[bgcolor=code]{c}
int 
f(int x, unsigned char y) {
   return (x+y)/2;
}
\end{minted}

\indent\text{Итак, опишем подходы чуть более подробно.}
\newline\indent\text{Подход с сохранением без знака требует приведения двух меньших типов без знака} \newline \text{к unsigned int. Это простое правило, и оно дает тип, который не зависит от среды} \newline \text{выполнения.}
\newline\indent\text{Подход с сохранением значений требует приведения малых типов к signed int,} \newline \text{если этот тип может исправно вмещать все значения исходного типа, иначе} \newline \text{приводится к unsigned int.} \newline \indent \text{Если среда выполнения представляет short собой нечто меньшее, чем int,} \newline \text{тогда unsigned short становится int; в противном случае становится unsigned int.} \newline \text{Обе схемы дают один и тот же ответ в подавляющем большинстве случаев, и обе} \newline \text{дают один и тот же эффективный результат в большинстве текущих реализаций.} \newline \text{В таких реализациях различия между ними проявляются только тогда, когда} \newline \text{выполняются следующие условия:
} \newline \text{ - Выражение, включающее unsigned char или unsigned short, дает int общий результат,} \newline \text{в котором установлен знаковый бит: т.е. либо унарная операция для такого типа,} \newline \text{либо двоичная операция, в которой другой операнд является int.

} \newline \text{ - Результат предыдущего выражения используется в контексте, в котором его} \newline \text{знаковость имеет значение:} \newline \text{
1. sizeof(int) < sizeof(long) и это ситуация, когда должно} \newline \text{произойти расширение до длинного типа} \newline \text{
2. это левый операнд оператора сдвига вправо (в реализации, где этот сдвиг} \newline \text{определяется как арифметический)} \newline \text{ 3. 
это либо операнды /, \%, <,<=,>, или >=.} \newline\indent \text{
В таких обстоятельствах возникает двусмысленность толкования.} \newline \text{Одним из важных результатов изучения этой проблемы является понимание того,} \newline \text{что высококачественным компиляторам может быть полезно искать такой} \newline \text{неоднозначно интерпретируемый код и предлагать (необязательно) диагностику.} \newline \indent\text{Правила сохранения без знака значительно увеличивают количество ситуаций,} \newline \text{которые приводят к неоднозначному результату, в то время как правила сохранения} \newline \text{значений сводят к минимуму такие ситуации. Таким образом, правила сохранения} \newline \text{значения считались более безопасными для начинающих или неосторожных} \newline \text{программистов. После долгих обсуждений Комитет принял решение в пользу} \newline \text{\underline{правил сохранения значения}, несмотря на то, что компиляторы UNIX} \newline \text{развивались в направлении сохранения без знака.}
\\
\newline\underline{Итак, ответ на поставленный вопрос задачи: сейчас все компиляторы расширяют} \newline \underline{малые целые следующим образом: делается попытка расширить до int, и только в} \newline \underline{случае неудачи - до unsigned int.} 
\newpage

\section*{Литература и веб-ресурсы}
\addcontentsline{toc}{section}{Литература и веб-ресурсы}
\begin{itemize}
\item[$-$]
\text{Н.В.Вдовикина, И.В.Машечкин, А.Н.Терехин, В.В.Тюляева} \newline \textbf{Программирование в ОС UNIX на языке Си}
\item[$-$]
\textbf{Rationale for
American National Standard for Information Systems -
      Programming Language -
C}
\item[$-$]
\text{microsin.net} \newline\textbf{Отличия кода ANSI C и кода K&R C}

\end{itemize}
\newpage
\end{document}
