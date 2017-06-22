#ifndef PARSER
#define PARSER

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <mem.h>
#include <string>
#include <string.h>
#include <math.h>
#include <set>
#include <map>
#include <sstream>
#include <cctype>
#include <list>
#include <vector>
#include <stack>
#include <deque>
#include <fstream>

//   -*-   -*-   -*-

enum E_TOKEN {
    E_OR, E_AND, E_EQ, E_NEQ, E_LT, E_LE, E_GE, E_GT, E_MINUS, E_PLUS, E_MUL, E_DIV, E_MOD, E_POW, E_NOT,
    E_UNARY_MINUS, E_UNARY_PLUS, E_FACTORIAL, E_LEFT_PAR, E_COMMA, E_RIGHT_PAR, E_NUM, E_VAR, E_ABS, E_ACOS, E_ACOSH,
    E_ACOT, E_ACOTH, E_ASIN, E_ASINH, E_ATAN, E_ATANH, E_BINOM, E_CEIL, E_CLAMP, E_COS, E_COSH, E_COT, E_COTH, E_E,
    E_EXP, E_FLOOR, E_IF, E_INF, E_ISFINITE, E_LN, E_LOG, E_MAX, E_MIN, E_PI, E_PROD, E_ROUND, E_SGN, E_SIN, E_SINH,
    E_SQRT, E_SUM, E_TAN, E_TANH, E_TRUNC,
    
    E_COLON, E_TOWARDS, E_ASSIGN
};

enum E_OPERATOR_ASSOCIATIVITY { E_LEFT, E_RIGHT };

enum E_TOKEN_TYPE { E_OPERATOR, E_CONTROL, E_NUMBER, E_VARIABLE, E_FUNCTION };

enum E_NUMBER_TYPE { E_REAL, E_INTEGER, E_BOOLEAN };

//   -*-   -*-   -*-

struct SOperatorInfo {
    int arity;
    int priority;
    E_OPERATOR_ASSOCIATIVITY associativity;
    std::string str;

    SOperatorInfo() {
    }

    SOperatorInfo(int arity, int priority, E_OPERATOR_ASSOCIATIVITY associativity, std::string str) :
        arity(arity), priority(priority), associativity(associativity), str(str)
    {
    }
};

struct STokenInfo {
    E_TOKEN token;
    E_TOKEN_TYPE type;
    std::string str;
    E_NUMBER_TYPE numType;

    STokenInfo() {
    }

    STokenInfo(E_TOKEN token, E_TOKEN_TYPE type, std::string str, E_NUMBER_TYPE numType) :
        token(token), type(type), str(str), numType(numType)
    {
    }
};

//   -*-   -*-   -*-

extern std::map<std::string, E_TOKEN> reserved;
extern std::map<E_TOKEN, SOperatorInfo> operators;

//   -*-   -*-   -*-

void BuildTokensMaps();
std::vector<STokenInfo> Tokenize(std::string formula);

//   -*-   -*-   -*-

extern int aVarInd;

//   -*-   -*-   -*-

int GetNextVar();
std::string MakeIndent(std::string str, char c, int size);
bool ParseFunction(std::vector<STokenInfo> &tokens, int &pos, std::string &code, int &varRes);
bool ParseExpression(std::vector<STokenInfo> &tokens, int &pos, std::string &code, int &varRes);

//   -*-   -*-   -*-

enum E_ARGUMENT_TYPE { E_ARGUMENT_TYPE_REAL, E_ARGUMENT_TYPE_REAL_VECTOR, E_ARGUMENT_TYPE_INTEGER, E_ARGUMENT_TYPE_INTEGER_VECTOR };

struct SArgumentInfo {
	E_ARGUMENT_TYPE type;
	int length; // for vectors
	bool varLength; // for vectors
	int firstVectorArgWithTheSameParamInd; // for vectors with variable length
};

struct SFunctionInfo {
	std::string name;
	std::vector<SArgumentInfo> args;
	int numberOfVectorArgs;
	SArgumentInfo retVal;
};

//   -*-   -*-   -*-

bool ParseBody(std::vector<STokenInfo> &tokens, int &pos, std::string &code, int &varRes);

//   -*-   -*-   -*-

#endif