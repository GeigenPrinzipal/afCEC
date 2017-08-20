\name{UserDefinedFormulas}
\alias{UserDefinedFormulas}

\title{
User-defined formulas
}

\section{Operators:}{
\subsection{List of the supported operators (from highest to the lowest priority)}{
\itemize{
\item{Arithmetic ! (factorial);}
\item{Arithmetic ^ (power);}
\item{Arithmetic unary +, - ;
\cr
Logical ! (not);
}
\item{Arithmetic *, /, \% (modulo division);}
\item{Arithmetic +, - ;}
\item{Relational <=, <, >, >= ;}
\item{Relational ==, != ;}
\item{Logical && (and);}
\item{Logical || (or);}
\item{Ternary conditional a ? b : c ;}
}
}

\subsection{Remarks}{
\itemize{
\item{! Computes the factorial of the argument of type integer. If the argument is not an integer value, it's rounded to the nearest one;}
\item{The relational operators produce the expression of the boolean type;}
\item{The logical operators can be used only with expressions of type boolean. Otherwise the exception is thrown and the \code{\link[=afCEC]{afCEC}} function returns with errors;}
\item{The boolean expressions can be used as arguments to the arithmetic operators. In that case they are implicitly converted to the 1 value (if the expression evaluates to TRUE) and 0 respectively, otherwise.
\cr
\cr
\bold{Example:}
\cr
\cr
The following expressions defined in three different ways are equivalent:
\cr
\cr
((x<=-1)*((x+1)^2+1))+((x>=1)*((x-1)^2+1))+(((x>-1)&&(x<1))*1)
\cr
\cr
x<=-1?(x+1)^2+1:x>=1?(x-1)^2+1:1
\cr
\cr
if(x<=-1,(x+1)^2+1,if(x>=1,(x-1)^2+1,1))
\cr
}
\item{The same applies for passing the boolean expressions as the function's real value arguments;}
}
}
}

\section{Built-in constants:}{
\itemize{
\item{e - Base of the natural logarithm;}
\item{pi - The pi value;}
}
}

\section{Built-in simple functions:}{
\itemize{
\item{abs(x) - Computes the absolute value of x;}
\item{acos(x) - Computes the inverse cosine of x and returns the value expressed in radians;}
\item{acosh(x) - Computes the inverse hyperbolic cosine of x and returns the value expressed in radians;}
\item{acot(x) - Computes the inverse cotangent of x and returns the value expressed in radians;}
\item{acoth(x) - Computes the inverse hyperbolic cotangent of x and returns the value expressed in radians;}
\item{asin(x) - Computes the inverse sine of x and returns the value expressed in radians;}
\item{asinh(x) - Computes the inverse hyperbolic sine of x and returns the value expressed in radians;}
\item{atan(x) - Computes the inverse tangent of x and returns the value expressed in radians;}
\item{atanh(x) - Computes the inverse hyperbolic tangent of x and returns the value expressed in radians;}
\item{binom(r,k) - Computes the r choose k, where r is an arbitrary real number and k is the integer. If k is not an integer value, it's rounded to the nearest one;}
\item{ceil(x) - Computes the ceiling function of x;}
\item{clamp(x,a,b) - Clamps x to the range [a,b];}
\item{cos(x) - Computes the cosine of x expressed in radians;}
\item{cosh(x) - Computes the hyperbolic cosine of x expressed in radians;}
\item{cot(x) - Computes the cotangent of tx expressed in radians;}
\item{coth(x) - Computes the hyperbolic cotangent of x expressed in radians;}
\item{exp(x) - Computes the exponent of x;}
\item{floor(x) - Computes the floor function of x;}
\item{frac(x) - Computes the fractional part of x;}
\item{isfinite(x) - Returns logical TRUE if x is the finite value. Otherwise, the logical FALSE is returned;}
\item{ln(x) - Computes the natural logarithm of x;}
\item{log(x) - Computes the 10-base logarithm of x;}
\item{log(a,x) - Computes the a-base logarithm of x;}
\item{round(x) - Rounds x to the nearest integer value;}
\item{sgn(x) - Computes the signum function of x;
\cr
\cr
If x < 0, sgn(x) returns -1;\cr
If x = 0, sgn(x) returns 0;\cr
If x > 0, sgn(x) returns 1;\cr
}
\item{sin(x) - Computes the sine of x expressed in radians;}
\item{sinh(x) - Computes the hyperbolic sine of x expressed in radians;}
\item{sqrt(x) - Computes the square root of x;}
\item{tan(x) - Computes the tangent of x expressed in radians;}
\item{tanh(x) - Computes the hyperbolic tangent of x expressed in radians;}
\item{trunc(x) - Truncates x (i.e., rounds x towards zero);}
}
}

\section{Built-in complex functions:}{
\itemize{

\item{if(<condition>, <expression_if_true>, <expression_if_false>) - Returns the expression determined by the <expression_if_true> argument, if the <condition> expression evaluates to TRUE. Otherwise, <expression_if_false> is returned. <condition> expression must have the boolean value. Otherwise the exception is thrown and the \code{\link[=afCEC]{afCEC}} function returns with errors.
\cr
\cr
\bold{Example:}
\cr
\cr
if(x<=-1,(x+1)^2+1,if(x>=1,(x-1)^2+1,1))
\cr
}

\item{max(x1, x2, ...) - Returns the maximum of the arguments. If no argument is specified (i.e., we have max()), the function returns the neutral element of the max operator (i.e., -INFINITY);}

\item{max(<index_variable>, <low_range>, <high_range>, <expression>) - Returns the maximum value of the expression determined by the <expression> argument depending on the index variable <index_variable> in the range [<low_range>, <high_range>]. <index_variable> is the indentifier of the index variable, while <low_range> and <high_range> are the integer values determining the range boundaries. If <low_range> or <high_range> are not integer values, they are rounded to the nearest ones. In the case where the interval [<low_range>, <high_range>] is degenerate (i.e., <high_range> < <low_range>), the function returns the neutral element of the max operator (i.e., -INFINITY);
\cr
\cr
\bold{Example:}
\cr
\cr
max(i,0,10,max(j,0,10,sin(i^j)))
\cr
}

\item{min(x1, x2, ...) - Returns the minimum of the arguments. If no argument is specified (i.e., we have min()), the function returns the neutral element of the min operator (i.e., +INFINITY);}

\item{min(<index_variable>, <low_range>, <high_range>, <expression>) - Returns the minimum value of the expression determined by the <expression> argument depending on the index variable <index_variable> in the range [<low_range>, <high_range>]. <index_variable> is the indentifier of the index variable, while <low_range> and <high_range> are the integer values determining the range boundaries. If <low_range> or <high_range> are not integer values, they are rounded to the nearest ones. In the case where the interval [<low_range>, <high_range>] is degenerate (i.e., <high_range> < <low_range>), the function returns the neutral element of the min operator (i.e., +INFINITY);
\cr
\cr
\bold{Example:}
\cr
\cr
min(i,0,10,min(j,0,10,sin(i^j)))
\cr
}

\item{prod(<index_variable>, <low_range>, <high_range>, <expression>) - Returns the product of the expression determined by the <expression> argument depending on the index variable <index_variable> over the range [<low_range>, <high_range>]. <index_variable> is the indentifier of the index variable, while <low_range> and <high_range> are the integer values determining the range boundaries. If <low_range> or <high_range> are not integer values, they are rounded to the nearest ones. In the case where the interval [<low_range>, <high_range>] is degenerate (i.e., <high_range> < <low_range>), the function returns the neutral element of the * operator (i.e., 1);
\cr
\cr
\bold{Example:}
\cr
\cr
prod(i,0,10,prod(j,0,10,sin(i^j)))
\cr
}

\item{sum(<index_variable>, <low_range>, <high_range>, <expression>) - Returns the sum of the expression determined by the <expression> argument depending on the index variable <index_variable> over the range [<low_range>, <high_range>]. <index_variable> is the indentifier of the index variable, while <low_range> and <high_range> are the integer values determining the range boundaries. If <low_range> or <high_range> are not integer values, they are rounded to the nearest ones. In the case where the interval [<low_range>, <high_range>] is degenerate (i.e., <high_range> < <low_range>), the function returns the neutral element of the + operator (i.e., 0);
\cr
\cr
\bold{Example:}
\cr
\cr
sum(i,0,10,sum(j,0,10,sin(i^j)))
}

}

}

\section{Declaring functions:}{
The function declaration consists of the function signature and the function formula placed after the signature.
\cr

\subsection{Function signature}{
\cr
\cr
The function signature is the expression of the form:
\cr
\cr
<funcion_name> : <domain_1> x ... x <domain_n> -> <counterdomain>,
\cr
\cr
where:
\cr
\cr
<function_name> - Is the identifier of the function;
\cr
\cr
<domain_i> - Specifies the domain of the i-th argument of the function. It can take the following values:
\itemize{
\item{R - Indicates that the function argument will be the real value scalar;}
\item{R^<integer_type_constant_expression> - Indicates that the function argument will be the real fixed size vector with number of elements equal to the integer value constant expression. If the expression is not of the type integer, it's rounded to the nearest integer value. The resulting number of elements must be the natural number > 0. Otherwise the exception is thrown and the \code{\link[=afCEC]{afCEC}} function returns with errors.}
\item{R^<identifier> - Indicates that the function argument will be the real vector with variable number of arguments accesible (inside the function formula) by the variable <indentifier>;}
\item{Z - Indicates that the function argument will be the integer value scalar;}
\item{Z^<integer_type_constant_expression> - Indicates that the function argument will be the integer fixed size vector with number of elements equal to the integer value constant expression. If the expression is not of the type integer, it's rounded to the nearest integer value. The resulting number of elements must be the natural number > 0. Otherwise the exception is thrown and the \code{\link[=afCEC]{afCEC}} function returns with errors.}
\item{Z^<identifier> - Indicates that the function argument will be the integer vector with variable number of arguments accesible (inside the function formula) by the variable <indentifier>;}
}
<counterdomain> - Specifies the counterdomain of the function. It can take the following values:
\itemize{
\item{R - Indicates that the counterdomain will be the real value scalar;}
\item{R^<integer_type_expression> - Indicates that the counterdomain will be the real vector with the number of elements given by the integer value expression depending on the variables indicating the size of the particular variable size vector arguments. If the expression is not of the type integer, it's rounded to the nearest integer value. The resulting number of elements must be the natural number > 0. Otherwise the exception is thrown and the \code{\link[=afCEC]{afCEC}} function returns with errors.
\cr
\cr
\bold{Example:}
\cr
\cr
f:R^n -> R^2*n+1
\cr
}
\item{Z - Indicates that the counterdomain will be the integer value scalar;}
\item{Z^<integer_type_expression> - Indicates that the counterdomain will be the integer vector with the number of elements given by the integer value expression depending on the variables indicating the size of the particular variable size vector arguments. If the expression is not of the type integer, it's rounded to the nearest integer value. The resulting number of elements must be the natural number > 0. Otherwise the exception is thrown and the \code{\link[=afCEC]{afCEC}} function returns with errors.
\cr
\cr
\bold{Example:}
\cr
\cr
f:R^n -> Z^2*n+1
\cr
}
}
\bold{Remarks:}
\itemize{
\item{The arguments passed to the function are strictly checked with the function signature for the size constraints, i.e., if the vector with more than one element is passed for the scalar argument (or vice-versa) or the vector with different number of elements than those declared in the signature in the case of the fixed size vector is passed as the fixed size vector argument, the exception is thrown and the \code{\link[=afCEC]{afCEC}} function returns with errors;}
\item{As it was implicitly stated above, one can pass the scalar argument if the fixed size vector with one element or the variable size vector is required. The same applies to passing the vectors containing exactly one element as the scalar arguments;}
\item{Conversions between the R and Z domain are allowed and performed automatically. In the case of the R to Z conversion, the rounding to the nearest integer value is used;}
\item{In the case of the variable size vectors, the same identifier may be used for more than one argument. In that case, the constraint on the size of that vectors is imposed in the sense that if there exist two or more vectors in the group with the same identifier that have different number of elements, then the exception is thrown and the \code{\link[=afCEC]{afCEC}} function returns with errors;
\cr
\cr
\bold{Example:}
\cr
\cr
f:R^n x R^m x R^n x R^m -> R
\cr
\cr
In that case, the first and third and the second and fourth vector must have exactly the same number of elements. Otherwise, the exception is thrown and the \code{\link[=afCEC]{afCEC}} function returns with errors;
}
\item{In the case of the fixed size vector arguments declarations using the integer type constant expressions, the globally visible constants passed to the \code{\link[=afCEC]{afCEC}} function via the \code{params} argument are also available to use in the expression. Similarly, the globally visible constants are also available in the integer type expressions in the declaration of the counterdomain;
\cr
\cr
\bold{Example:}
\cr
\cr
Suppose, that the globally visible constant d indicating the degree of the multivariate polynomial was passed to the function \code{\link[=afCEC]{afCEC}}. The declaration of that polynomial may look like as follows:
\cr
\cr
f:R^n -> R^n^d
}
}
}



}
