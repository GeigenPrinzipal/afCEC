#include ".\Parser.h"

//   -*-   -*-   -*-

std::map<std::string, E_TOKEN> reserved;
std::map<E_TOKEN, SOperatorInfo> operators;
int aVarInd = 0;

//   -*-   -*-   -*-

void BuildTokensMaps() {
    std::map<std::string, E_TOKEN>::iterator it1 = reserved.begin();
    it1 = reserved.insert(it1, std::pair<std::string, E_TOKEN>("abs", E_ABS));
    it1 = reserved.insert(it1, std::pair<std::string, E_TOKEN>("acos", E_ACOS));
    it1 = reserved.insert(it1, std::pair<std::string, E_TOKEN>("acosh", E_ACOSH));
    it1 = reserved.insert(it1, std::pair<std::string, E_TOKEN>("acot", E_ACOT));
    it1 = reserved.insert(it1, std::pair<std::string, E_TOKEN>("acoth", E_ACOTH));
    it1 = reserved.insert(it1, std::pair<std::string, E_TOKEN>("asin", E_ASIN));
    it1 = reserved.insert(it1, std::pair<std::string, E_TOKEN>("asinh", E_ASINH));
    it1 = reserved.insert(it1, std::pair<std::string, E_TOKEN>("atan", E_ATAN));
    it1 = reserved.insert(it1, std::pair<std::string, E_TOKEN>("atanh", E_ATANH));
    it1 = reserved.insert(it1, std::pair<std::string, E_TOKEN>("binom", E_BINOM));
    it1 = reserved.insert(it1, std::pair<std::string, E_TOKEN>("ceil", E_CEIL));
    it1 = reserved.insert(it1, std::pair<std::string, E_TOKEN>("clamp", E_CLAMP));
    it1 = reserved.insert(it1, std::pair<std::string, E_TOKEN>("cos", E_COS));
    it1 = reserved.insert(it1, std::pair<std::string, E_TOKEN>("cosh", E_COSH));
    it1 = reserved.insert(it1, std::pair<std::string, E_TOKEN>("cot", E_COT));
    it1 = reserved.insert(it1, std::pair<std::string, E_TOKEN>("coth", E_COTH));
    it1 = reserved.insert(it1, std::pair<std::string, E_TOKEN>("e", E_E));
    it1 = reserved.insert(it1, std::pair<std::string, E_TOKEN>("exp", E_EXP));
    it1 = reserved.insert(it1, std::pair<std::string, E_TOKEN>("floor", E_FLOOR));
    it1 = reserved.insert(it1, std::pair<std::string, E_TOKEN>("if", E_IF));
    it1 = reserved.insert(it1, std::pair<std::string, E_TOKEN>("inf", E_INF));
    it1 = reserved.insert(it1, std::pair<std::string, E_TOKEN>("isfinite", E_ISFINITE));
    it1 = reserved.insert(it1, std::pair<std::string, E_TOKEN>("ln", E_LN));
    it1 = reserved.insert(it1, std::pair<std::string, E_TOKEN>("log", E_LOG));
    it1 = reserved.insert(it1, std::pair<std::string, E_TOKEN>("max", E_MAX));
    it1 = reserved.insert(it1, std::pair<std::string, E_TOKEN>("min", E_MIN));
    it1 = reserved.insert(it1, std::pair<std::string, E_TOKEN>("pi", E_PI));
    it1 = reserved.insert(it1, std::pair<std::string, E_TOKEN>("prod", E_PROD));
    it1 = reserved.insert(it1, std::pair<std::string, E_TOKEN>("round", E_ROUND));
    it1 = reserved.insert(it1, std::pair<std::string, E_TOKEN>("sgn", E_SGN));
    it1 = reserved.insert(it1, std::pair<std::string, E_TOKEN>("sin", E_SIN));
    it1 = reserved.insert(it1, std::pair<std::string, E_TOKEN>("sinh", E_SINH));
    it1 = reserved.insert(it1, std::pair<std::string, E_TOKEN>("sqrt", E_SQRT));
    it1 = reserved.insert(it1, std::pair<std::string, E_TOKEN>("sum", E_SUM));
    it1 = reserved.insert(it1, std::pair<std::string, E_TOKEN>("tan", E_TAN));
    it1 = reserved.insert(it1, std::pair<std::string, E_TOKEN>("tanh", E_TANH));
    reserved.insert(it1, std::pair<std::string, E_TOKEN>("trunc", E_TRUNC));
    
    std::map<E_TOKEN, SOperatorInfo>::iterator it2 = operators.begin();
    it2 = operators.insert(it2, std::pair<E_TOKEN, SOperatorInfo>(E_OR, SOperatorInfo(2, 0, E_LEFT, "||")));
    it2 = operators.insert(it2, std::pair<E_TOKEN, SOperatorInfo>(E_AND, SOperatorInfo(2, 1, E_LEFT, "&&")));
    it2 = operators.insert(it2, std::pair<E_TOKEN, SOperatorInfo>(E_EQ, SOperatorInfo(2, 2, E_LEFT, "==")));
    it2 = operators.insert(it2, std::pair<E_TOKEN, SOperatorInfo>(E_NEQ, SOperatorInfo(2, 2, E_LEFT, "!=")));
    it2 = operators.insert(it2, std::pair<E_TOKEN, SOperatorInfo>(E_LT, SOperatorInfo(2, 3, E_LEFT, "<")));
    it2 = operators.insert(it2, std::pair<E_TOKEN, SOperatorInfo>(E_LE, SOperatorInfo(2, 3, E_LEFT, "<=")));
    it2 = operators.insert(it2, std::pair<E_TOKEN, SOperatorInfo>(E_GE, SOperatorInfo(2, 3, E_LEFT, ">=")));
    it2 = operators.insert(it2, std::pair<E_TOKEN, SOperatorInfo>(E_GT, SOperatorInfo(2, 3, E_LEFT, ">")));
    it2 = operators.insert(it2, std::pair<E_TOKEN, SOperatorInfo>(E_MINUS, SOperatorInfo(2, 4, E_LEFT, "-")));
    it2 = operators.insert(it2, std::pair<E_TOKEN, SOperatorInfo>(E_PLUS, SOperatorInfo(2, 4, E_LEFT, "+")));
    it2 = operators.insert(it2, std::pair<E_TOKEN, SOperatorInfo>(E_MUL, SOperatorInfo(2, 5, E_LEFT, "*")));
    it2 = operators.insert(it2, std::pair<E_TOKEN, SOperatorInfo>(E_DIV, SOperatorInfo(2, 5, E_LEFT, "/")));
    it2 = operators.insert(it2, std::pair<E_TOKEN, SOperatorInfo>(E_MOD, SOperatorInfo(2, 5, E_LEFT, "%")));
    it2 = operators.insert(it2, std::pair<E_TOKEN, SOperatorInfo>(E_POW, SOperatorInfo(2, 6, E_RIGHT, "^")));
    it2 = operators.insert(it2, std::pair<E_TOKEN, SOperatorInfo>(E_NOT, SOperatorInfo(1, 7, E_RIGHT, "!")));
    it2 = operators.insert(it2, std::pair<E_TOKEN, SOperatorInfo>(E_UNARY_MINUS, SOperatorInfo(1, 7, E_RIGHT, "-")));
    it2 = operators.insert(it2, std::pair<E_TOKEN, SOperatorInfo>(E_UNARY_PLUS, SOperatorInfo(1, 7, E_RIGHT, "+")));
    operators.insert(it2, std::pair<E_TOKEN, SOperatorInfo>(E_FACTORIAL, SOperatorInfo(1, 8, E_LEFT, "!")));
}
   
//   -*-   -*-   -*-

std::vector<STokenInfo> Tokenize(std::string formula) {
    std::istringstream iss(formula);
    std::list<STokenInfo> lst;
    bool binary = false;
    while (!iss.eof()) {
        std::string str;
        iss >> str;
        int j = 0;
        while (j < str.length()) {
            switch (str[j]) {
                case ':' : {
                    std::string substr = str.substr(j, 1);
                    lst.push_back(STokenInfo(E_COLON, E_CONTROL, substr, ((E_NUMBER_TYPE)0)));
                    ++j;
                    break;
                }
                case '|' : {
                    if ((j + 1 < str.length()) && (str[j + 1] == '|')) {
                        std::string substr = str.substr(j, 2);
                        lst.push_back(STokenInfo(E_OR, E_OPERATOR, substr, ((E_NUMBER_TYPE)0)));
                        binary = false;
                        j += 2;
                    } else
                        throw 0;
                    break;
                }
                case '&' : {
                    if ((j + 1 < str.length()) && (str[j + 1] == '&')) {
                        std::string substr = str.substr(j, 2);
                        lst.push_back(STokenInfo(E_AND, E_OPERATOR, substr, ((E_NUMBER_TYPE)0)));
                        binary = false;
                        j += 2;
                    } else
                        throw 0;
                    break;
                }
                case '=' : {
                    if ((j + 1 < str.length()) && (str[j + 1] == '=')) {
                        std::string substr = str.substr(j, 2);
                        lst.push_back(STokenInfo(E_EQ, E_OPERATOR, substr, ((E_NUMBER_TYPE)0)));
                        binary = false;
                        j += 2;
                    } else {
                        std::string substr = str.substr(j, 1);
                        lst.push_back(STokenInfo(E_ASSIGN, E_CONTROL, substr, ((E_NUMBER_TYPE)0)));
                        ++j;
                    }
                    break;
                }
                case '!' : {
                    if ((j + 1 < str.length()) && (str[j + 1] == '=')) {
                        std::string substr = str.substr(j, 2);
                        lst.push_back(STokenInfo(E_NEQ, E_OPERATOR, substr, ((E_NUMBER_TYPE)0)));
                        binary = false;
                        j += 2;
                    } else {
                        std::string substr = str.substr(j, 1);
                        if (binary) lst.push_back(STokenInfo(E_FACTORIAL, E_OPERATOR, substr, ((E_NUMBER_TYPE)0)));
                        else {
                            lst.push_back(STokenInfo(E_NOT, E_OPERATOR, substr, ((E_NUMBER_TYPE)0)));
                            binary = false;
                        }
                        ++j;
                    }
                    break;
                }
                case '<' : {
                    if ((j + 1 < str.length()) && (str[j + 1] == '=')) {
                        std::string substr = str.substr(j, 2);
                        lst.push_back(STokenInfo(E_LE, E_OPERATOR, substr, ((E_NUMBER_TYPE)0)));
                        j += 2;
                    } else {
                        std::string substr = str.substr(j, 1);
                        lst.push_back(STokenInfo(E_LT, E_OPERATOR, substr, ((E_NUMBER_TYPE)0)));
                        ++j;
                    }
                    binary = false;
                    break;
                }
                case '>' : {
                    if ((j + 1 < str.length()) && (str[j + 1] == '=')) {
                        std::string substr = str.substr(j, 2);
                        lst.push_back(STokenInfo(E_GE, E_OPERATOR, substr, ((E_NUMBER_TYPE)0)));
                        j += 2;
                    } else {
                        std::string substr = str.substr(j, 1);
                        lst.push_back(STokenInfo(E_GT, E_OPERATOR, substr, ((E_NUMBER_TYPE)0)));
                        ++j;
                    }
                    binary = false;
                    break;
                }
                case '-' : {
                    if ((j + 1 < str.length()) && (str[j + 1] == '>')) {
                        std::string substr = str.substr(j, 2);
                        lst.push_back(STokenInfo(E_TOWARDS, E_CONTROL, substr, ((E_NUMBER_TYPE)0)));
                        j += 2;
                    } else {
                        std::string substr = str.substr(j, 1);
                        if (binary) {
                            lst.push_back(STokenInfo(E_MINUS, E_OPERATOR, substr, ((E_NUMBER_TYPE)0)));
                            binary = false;
                        } else
                            lst.push_back(STokenInfo(E_UNARY_MINUS, E_OPERATOR, substr, ((E_NUMBER_TYPE)0)));
                        ++j;
                    }
                    break;
                }
                case '+' : {
                    std::string substr = str.substr(j, 1);
                    if (binary) {
                        lst.push_back(STokenInfo(E_PLUS, E_OPERATOR, substr, ((E_NUMBER_TYPE)0)));
                        binary = false;
                    } else
                        lst.push_back(STokenInfo(E_UNARY_PLUS, E_OPERATOR, substr, ((E_NUMBER_TYPE)0)));
                    ++j;
                    break;
                }
                case '*' : {
                    std::string substr = str.substr(j, 1);
                    lst.push_back(STokenInfo(E_MUL, E_OPERATOR, substr, ((E_NUMBER_TYPE)0)));
                    binary = false;
                    ++j;
                    break;
                }
                case '/' : {
                    std::string substr = str.substr(j, 1);
                    lst.push_back(STokenInfo(E_DIV, E_OPERATOR, substr, ((E_NUMBER_TYPE)0)));
                    binary = false;
                    ++j;
                    break;
                }
                case '%' : {
                    std::string substr = str.substr(j, 1);
                    lst.push_back(STokenInfo(E_MOD, E_OPERATOR, substr, ((E_NUMBER_TYPE)0)));
                    binary = false;
                    ++j;
                    break;
                }
                case '^' : {
                    std::string substr = str.substr(j, 1);
                    lst.push_back(STokenInfo(E_POW, E_OPERATOR, substr, ((E_NUMBER_TYPE)0)));
                    binary = false;
                    ++j;
                    break;
                }
                case '(' : {
                    std::string substr = str.substr(j, 1);
                    lst.push_back(STokenInfo(E_LEFT_PAR, E_CONTROL, substr, ((E_NUMBER_TYPE)0)));
                    binary = false;
                    ++j;
                    break;
                }
                case ')' : {
                    std::string substr = str.substr(j, 1);
                    lst.push_back(STokenInfo(E_RIGHT_PAR, E_CONTROL, substr, ((E_NUMBER_TYPE)0)));
                    binary = true;
                    ++j;
                    break;
                }
                case ',' : {
                    std::string substr = str.substr(j, 1);
                    lst.push_back(STokenInfo(E_COMMA, E_CONTROL, substr, ((E_NUMBER_TYPE)0)));
                    binary = false;
                    ++j;
                    break;
                }
                default : {
                    if ((isdigit(str[j]) != 0) || (str[j] == '.')) {
                        int k = j + 1;
                        if (str[j] != '.') {
                            while ((k < str.length()) && (isdigit(str[k]) != 0)) ++k;
                            if ((k < str.length()) && (str[k] == '.')) ++k;
                            while ((k < str.length()) && (isdigit(str[k]) != 0)) ++k;
                        } else {
                            if ((k < str.length()) && (isdigit(str[k]) != 0)) {
                                ++k;
                                while ((k < str.length()) && (isdigit(str[k]) != 0)) ++k;
                            } else
                                throw 0;
                        }
                        int l = k;
                        if ((l < str.length()) && ((str[l] == 'e') || (str[l] == 'E'))) {
                            ++l;
                            if ((l < str.length()) && ((str[l] == '-') || (str[l] == '+'))) ++l;
                            if ((l < str.length()) && (isdigit(str[l]) != 0)) {
                                ++l;
                                while ((l < str.length()) && (isdigit(str[l]) != 0)) ++l;
                                k = l;
                            }
                        }
                        std::string substr = str.substr(j, k - j);
                        std::stringstream ss;
                        ss << substr;
                        double dVal;
                        ss >> dVal;
                        if (ss.fail()) throw 0;
                        ss.str("");
                        ss.clear();
                        int iVal = dVal;
                        ss << iVal;
                        if (ss.str() != substr) lst.push_back(STokenInfo(E_NUM, E_NUMBER, substr, E_REAL));
                        else {
                            if ((iVal != 0) && (iVal != 1)) lst.push_back(STokenInfo(E_NUM, E_NUMBER, substr, E_INTEGER));
                            else
                                lst.push_back(STokenInfo(E_NUM, E_NUMBER, substr, E_BOOLEAN));
                        }
                        binary = true;
                        j = k;
                    } else if ((isalpha(str[j]) != 0) || (str[j] == '_')) {
                        // [a-zA-Z_][a-zA-Z0-9_]* - regex for identifier
                        int k = j + 1;
                        while ((k < str.length()) && ((isalnum(str[k]) != 0) || (str[k] == '_'))) ++k;
                        std::string substr = str.substr(j, k - j);
                        std::map<std::string, E_TOKEN>::iterator it = reserved.find(substr);
                        if (it != reserved.end()) {
                            if ((substr == "e") || (substr == "pi") || (substr == "inf"))
                                lst.push_back(STokenInfo(E_NUM, E_NUMBER, substr, E_REAL));
                            else
                                lst.push_back(STokenInfo(it->second, E_FUNCTION, substr, ((E_NUMBER_TYPE)0)));
                        } else
                            lst.push_back(STokenInfo(E_VAR, E_VARIABLE, substr, ((E_NUMBER_TYPE)0)));
                        binary = true;
                        j = k;
                    } else
                        throw 0;
                }
            }
        }
    }
    std::vector<STokenInfo> tokens(lst.begin(), lst.end());
    return tokens;
}

//   -*-   -*-   -*-

int GetNextVar() {
    return aVarInd++;
}

//   -*-   -*-   -*-

STokenInfo GetNextToken(std::vector<STokenInfo> tokens, int &pos) {
    if (pos < tokens.size()) return tokens[pos];
    else
        throw 0;
}

//   -*-   -*-   -*-

std::string MakeIndent(std::string str, char c, int size) {
    std::string indent(size, c);
    std::string res = "";
    int lo = 0;
    for (int i = 0; i < str.length(); ++i) {
        if (str[i] == '\n') {
            res += indent + str.substr(lo, (i - lo) + 1);
            lo = i + 1;
        }
    }
    if (lo < str.length()) res += indent + str.substr(lo, str.length() - lo);
    return res;
}

//   -*-   -*-   -*-

bool ParseExpression(std::vector<STokenInfo> &tokens, int &pos, std::string &code, int &varRes) {
    try {
        std::stack<E_TOKEN> stOp;
        std::stack<int> stVar;
        std::ostringstream oss;
        bool done = false;
        
        while ((pos < tokens.size()) && (!done)) {
            STokenInfo token = GetNextToken(tokens, pos);
            switch (token.type) {
                case E_OPERATOR : {
                    if (!stOp.empty()) {
                        SOperatorInfo op = operators[token.token];
                        SOperatorInfo opSt = operators[stOp.top()];
                        while (
                            (!stOp.empty()) 
                            &&
                                (((opSt.associativity == E_LEFT) && (opSt.priority >= op.priority)) 
                                ||
                                    ((opSt.associativity == E_RIGHT) 
                                    && 
                                    (opSt.priority > op.priority) 
                                    && 
                                    (((stOp.top() != E_UNARY_MINUS) && (stOp.top() != E_UNARY_PLUS)) || (token.token != E_POW)))
                                )
                        ) {
                            if (opSt.arity == 1) {
                                int varRes1 = stVar.top();
                                stVar.pop();
                                int varNew1 = GetNextVar();
                                if (stOp.top() != E_FACTORIAL) oss << "double $" << varNew1 << " = " << opSt.str << "$" << varRes1 << ";\n";
                                else
                                    oss << "double $" << varNew1 << " = Factorial($" << varRes1 << ");\n";
                                stVar.push(varNew1);                     
                            } else {
                                int varRes2 = stVar.top();
                                stVar.pop();
                                int varRes1 = stVar.top();
                                stVar.top();
                                int varNew1 = GetNextVar();
                                if (stOp.top() != E_POW) oss << "double $" << varNew1 << " = $" << varRes1 << " " << opSt.str << " $" << varRes2 << ";\n";
                                else
                                    oss << "double $" << varNew1 << " = pow($" << varRes1 << ", $" << varRes2 << ");\n";
                                stVar.push(varNew1);
                            }
                            stOp.pop();
                            if (!stOp.empty()) opSt = operators[stOp.top()];
                        }
                    }
                    stOp.push(token.token);
                    ++pos;
                    break;
                }
                case E_CONTROL : {
                    switch (token.token) {
                        case E_LEFT_PAR : {
                            std::string code1;
                            int varRes1;
                            
                            ++pos;
                            if (!ParseExpression(tokens, pos, code1, varRes1)) throw 0;
                            if (GetNextToken(tokens, pos).token != E_RIGHT_PAR) throw 0;
                            ++pos;
                            if ((pos < tokens.size()) && (GetNextToken(tokens, pos).type != E_OPERATOR)) done = true;
                            
                            oss << code1;
                            stVar.push(varRes1);
                            break;
                        }
                        default : {
                            while (!stOp.empty()) {
                                E_TOKEN token = stOp.top();
                                stOp.pop();
                                SOperatorInfo opSt = operators[token];
                                if (opSt.arity == 1) {
                                    int varRes1 = stVar.top();
                                    stVar.pop();
                                    int varNew1 = GetNextVar();
                                    if (token != E_FACTORIAL) oss << "double $" << varNew1 << " = " << opSt.str << "$" << varRes1 << ";\n";
                                    else
                                        oss << "double $" << varNew1 << " = Factorial($" << varRes1 << ");\n";
                                    stVar.push(varNew1);
                                } else {                                    
                                    int varRes2 = stVar.top();
                                    stVar.pop();
                                    int varRes1 = stVar.top();
                                    stVar.pop();
                                    int varNew1 = GetNextVar();
                                    if (token != E_POW) oss << "double $" << varNew1 << " = $" << varRes1 << " " << opSt.str << " $" << varRes2 << ";\n";
                                    else
                                        oss << "double $" << varNew1 << " = pow($" << varRes1 << ", $" << varRes2 << ");\n";
                                    stVar.push(varNew1);
                                }
                            }
                            done = true;
                        }
                    }
                    break;
                }
                case E_NUMBER : {
                    ++pos;
                    int varNew1 = GetNextVar();
                    oss << "double $" << varNew1 << " = " << token.str << ";\n";
                    stVar.push(varNew1);
                    if ((pos < tokens.size()) && (GetNextToken(tokens, pos).type != E_OPERATOR)) done = true;
                    break;
                }
                case E_VARIABLE : {
                    /*
                        <expr>
                        $1 = <var>(<expr>.res);
                        
                        -----
                        
                        $1 = <var>;
                    */
                    
                    std::string code1;
                    int varRes1;
                    
                    ++pos;
                    if ((pos < tokens.size()) && (GetNextToken(tokens, pos).token == E_LEFT_PAR)) {
                        ++pos;
                        if (!ParseExpression(tokens, pos, code1, varRes1)) throw 0;
                        if (GetNextToken(tokens, pos).token != E_RIGHT_PAR) throw 0;
                        ++pos;
                        
                        int varNew1 = GetNextVar();
                        // <expr>
                        oss << code1;
                        // $1 = <var>(<expr>.res);
                        oss << "double $" << varNew1 << " = " << /*token.str*/ "$0" << "($" << varRes1 << " - 1);\n";
                        
                        stVar.push(varNew1);
                    } else {
                        int varNew1 = GetNextVar();
                        oss << "double $" << varNew1 << " = " << /*token.str*/ "$0" << ";\n";
                        stVar.push(varNew1);
                    }                    
                    
                    if ((pos < tokens.size()) && (GetNextToken(tokens, pos).type != E_OPERATOR)) done = true;
                    break;
                }
                case E_FUNCTION : {
                    std::string code1;
                    int varRes1;
                    
                    if (!ParseFunction(tokens, pos, code1, varRes1)) throw 0;
                    if ((pos < tokens.size()) && (GetNextToken(tokens, pos).type != E_OPERATOR)) done = true;
                    
                    oss << code1;
                    stVar.push(varRes1);
                    break;
                } 
            }
        }
        
        while (!stOp.empty()) {
            E_TOKEN token = stOp.top();
            stOp.pop();
            SOperatorInfo opSt = operators[token];
            if (opSt.arity == 1) {
                int varRes1 = stVar.top();
                stVar.pop();
                int varNew1 = GetNextVar();
                if (token != E_FACTORIAL) oss << "double $" << varNew1 << " = " << opSt.str << "$" << varRes1 << ";\n";
                else
                    oss << "double $" << varNew1 << " = Factorial($" << varRes1 << ");\n";
                stVar.push(varNew1);
            } else {                                    
                int varRes2 = stVar.top();
                stVar.pop();
                int varRes1 = stVar.top();
                stVar.pop();
                int varNew1 = GetNextVar();
                if (token != E_POW) oss << "double $" << varNew1 << " = $" << varRes1 << " " << opSt.str << " $" << varRes2 << ";\n";
                else
                    oss << "double $" << varNew1 << " = pow($" << varRes1 << ", $" << varRes2 << ");\n";
                stVar.push(varNew1);
            }
        }
        
        code = oss.str();
        varRes = stVar.top();
        stVar.pop();
    } catch (...) {
        return false;
    }
    return true;
}

//   -*-   -*-   -*-

bool ParseFunction(std::vector<STokenInfo> &tokens, int &pos, std::string &code, int &varRes) {
    int tabSize = 4; // !!!
    try {
        std::deque<int> dVar;
        std::ostringstream oss;
        
        STokenInfo token1 = GetNextToken(tokens, pos);
        ++pos;
        if (GetNextToken(tokens, pos).token != E_LEFT_PAR) throw 0;
        ++pos;
        switch (token1.token) {
            case E_IF : {
                break;
            }
            case E_MAX : {
                /*
                    max(<cv>, <lo>, <hi>, <expr>)
                    
                    <lo>
                    <hi>
                    double $1 = -datum::inf;
                    int <cv> = <lo>.res;
                    while (<cv> <= <hi>.res) {
                        <expr>
                        if (<expr>.res > $1) $1 = <expr>.res;
                        ++cv;
                    }
                */
                STokenInfo token2 = GetNextToken(tokens, pos);
                ++pos;
                
                std::string code1, code2, code3;
                int varRes1, varRes2, varRes3;
                
                if (GetNextToken(tokens, pos).token != E_COMMA) throw 0;
                ++pos;
                if (!ParseExpression(tokens, pos, code1, varRes1)) throw 0; // <lo>
                if (GetNextToken(tokens, pos).token != E_COMMA) throw 0;
                ++pos;
                if (!ParseExpression(tokens, pos, code2, varRes2)) throw 0; // <hi>
                if (GetNextToken(tokens, pos).token != E_COMMA) throw 0;
                ++pos;
                if (!ParseExpression(tokens, pos, code3, varRes3)) throw 0; // <expr>
                if (GetNextToken(tokens, pos).token != E_RIGHT_PAR) throw 0;
                ++pos;
                
                // <lo>
                oss << code1;
                // <hi>
                oss << code2;
                int varNew1 = GetNextVar();
                // double $1 = -datum::inf;
                oss << "double $" << varNew1 << " = -datum::inf;\n";
                // int <cv> = <lo>.res;
                oss << "int " << token2.str << " = $" << varRes1 << ";\n";
                // while (<cv> <= <hi>.res) {
                oss << "while (" << token2.str << " <= $" << varRes2 << ") {\n";
                // <expr>
                oss << MakeIndent(code3, ' ', tabSize);
                // if (<expr>.res > $1) $1 = <expr>.res;
                oss << std::string(tabSize, ' ') << "if ($" << varRes3 << " > $" << varNew1 << ") $" << varNew1 << " = $" << varRes3 << ";\n";
                // ++cv;
                oss << std::string(tabSize, ' ') << "++" << token2.str << ";\n";
                // }
                oss << "}\n";
                
                varRes = varNew1;
                break;
            }
            case E_MIN : {
                /*
                    max(<cv>, <lo>, <hi>, <expr>)
                    
                    <lo>
                    <hi>
                    double $1 = -datum::inf;
                    int <cv> = <lo>.res;
                    while (<cv> <= <hi>.res) {
                        <expr>
                        if (<expr>.res > $1) $1 = <expr>.res;
                        ++cv;
                    }
                */
                STokenInfo token2 = GetNextToken(tokens, pos);
                ++pos;
                
                std::string code1, code2, code3;
                int varRes1, varRes2, varRes3;
                
                if (GetNextToken(tokens, pos).token != E_COMMA) throw 0;
                ++pos;
                if (!ParseExpression(tokens, pos, code1, varRes1)) throw 0; // <lo>
                if (GetNextToken(tokens, pos).token != E_COMMA) throw 0;
                ++pos;
                if (!ParseExpression(tokens, pos, code2, varRes2)) throw 0; // <hi>
                if (GetNextToken(tokens, pos).token != E_COMMA) throw 0;
                ++pos;
                if (!ParseExpression(tokens, pos, code3, varRes3)) throw 0; // <expr>
                if (GetNextToken(tokens, pos).token != E_RIGHT_PAR) throw 0;
                ++pos;
                
                // <lo>
                oss << code1;
                // <hi>
                oss << code2;
                int varNew1 = GetNextVar();
                // double $1 = -datum::inf;
                oss << "double $" << varNew1 << " = datum::inf;\n";
                // int <cv> = <lo>.res;
                oss << "int " << token2.str << " = $" << varRes1 << ";\n";
                // while (<cv> <= <hi>.res) {
                oss << "while (" << token2.str << " <= $" << varRes2 << ") {\n";
                // <expr>
                oss << MakeIndent(code3, ' ', tabSize);
                // if (<expr>.res > $1) $1 = <expr>.res;
                oss << std::string(tabSize, ' ') << "if ($" << varRes3 << " < $" << varNew1 << ") $" << varNew1 << " = $" << varRes3 << ";\n";
                // ++cv;
                oss << std::string(tabSize, ' ') << "++" << token2.str << ";\n";
                // }
                oss << "}\n";
                
                varRes = varNew1;
                break;
            }
            case E_PROD : {
                break;
            }
            case E_SUM : {
                /*
                    sum(<cv>, <lo>, <hi>, <expr>)
                    
                    double $1;
                    {
                        <lo>
                        <hi>
                        $1 = 0;
                        int <cv> = $<lo>;
                        while (<cv> <= $<hi>) {
                            <expr>
                            $1 += $<expr>;
                            ++cv;
                        }
                    }
                */
                STokenInfo token = GetNextToken(tokens, pos);
                ++pos;
                
                std::string code1, code2, code3;
                int varRes1, varRes2, varRes3;
                
                if (GetNextToken(tokens, pos).token != E_COMMA) throw 0;
                ++pos;
                if (!ParseExpression(tokens, pos, code1, varRes1)) throw 0; // <lo>
                if (GetNextToken(tokens, pos).token != E_COMMA) throw 0;
                ++pos;
                if (!ParseExpression(tokens, pos, code2, varRes2)) throw 0; // <hi>
                if (GetNextToken(tokens, pos).token != E_COMMA) throw 0;
                ++pos;
                if (!ParseExpression(tokens, pos, code3, varRes3)) throw 0; // <expr>
                if (GetNextToken(tokens, pos).token != E_RIGHT_PAR) throw 0;
                ++pos;
                
                //double $1;
                int varNew = GetNextVar();
                oss << "double $" << varNew << ";" << std::endl;
                // {
                oss << "{" << std::endl;
                // <lo>
                oss << MakeIndent(code1, ' ', tabSize);
                // <hi>
                oss << MakeIndent(code2, ' ', tabSize);
                // $1 = 0;
                oss << std::string(tabSize, ' ') << "$" << varNew << " = 0;" << std::endl;
                // int <cv> = $<lo>;
                oss << std::string(tabSize, ' ') << "int " << token.str << " = $" << varRes1 << ";" << std::endl;
                // while (<cv> <= $<hi>) {
                oss << std::string(tabSize, ' ') << "while (" << token.str << " <= $" << varRes2 << ") {" << std::endl;
                // <expr>
                oss << MakeIndent(code3, ' ', tabSize * 2);
                // $1 += $<expr>;
                oss << std::string(tabSize * 2, ' ') << "$" << varNew << " += $" << varRes3 << ";" << std::endl;
                // ++cv;
                oss << std::string(tabSize * 2, ' ') << "++" << token.str << ";" << std::endl;
                // }
                oss << std::string(tabSize, ' ') << "}" << std::endl;
                // }
                oss << "}" << std::endl;
                
                varRes = varNew;
                break;
            }
            default : {
                bool done = false;
                int argsNum = 0;
                while (!done) {
                    std::string code1;
                    int varRes1;
                    
                    if (!ParseExpression(tokens, pos, code1, varRes1)) throw 0;
                    ++argsNum;
                    
                    oss << code1;
                    dVar.push_back(varRes1);
                    
                    STokenInfo token2 = GetNextToken(tokens, pos);
                    ++pos;
                    switch (token2.token) {
                        case E_RIGHT_PAR : {
                            done = true;
                            break;
                        }
                        case E_COMMA : {
                            break;
                        }
                        default : {
                            throw 0;
                        }
                    }
                }
                int varNew1 = GetNextVar();
                oss << "double $" << varNew1 << " = " << token1.str << "(";
                for (int i = 0; i < argsNum; ++i) {
                    oss << "$" << dVar.front();
                    dVar.pop_front();
                    if (i < argsNum - 1) oss << ", ";
                }
                oss << ");\n";
                varRes = varNew1;
            }
        }
        
        code = oss.str();
    } catch (...) {
        return false;
    }
    return true;
}

//   -*-   -*-   -*-















//   -*-   -*-   -*-

bool ParseBody(std::vector<STokenInfo> &tokens, int &pos, std::string &code, int &varRes) {
    const int tabSize = 4;
    try {
        bool done = false;
		std::ostringstream oss;
        
        while ((pos < tokens.size()) && (!done)) {
            std::map<std::string, int> symbols;
		    std::list<SArgumentInfo> args;
		    SFunctionInfo func;
            func.numberOfVectorArgs = 0;
            
            // Parsing the function domain and counterdomain specification of the form:
            // f : <dom1> x <dom2> x ... <dom_n> -> <dom_(n+1)>
            STokenInfo token = GetNextToken(tokens, pos);
            if (token.type == E_VARIABLE) func.name = token.str;
			else
				throw 0;
			++pos;
			if (GetNextToken(tokens, pos).token != E_COLON) throw 0;
			++pos;
			bool done = false;
			
			do {
				STokenInfo token = GetNextToken(tokens, pos);
				++pos;
                SArgumentInfo arg;
				switch (token.str[0]) {
					case 'R' : {
						if (GetNextToken(tokens, pos).str[0] == '^') {
							++pos;
							STokenInfo token = GetNextToken(tokens, pos);
							++pos;
							switch (token.type) {
                                case E_VARIABLE : {
                                    arg.varLength = true;
                                    if (symbols.find(token.str) != symbols.end()) arg.firstVectorArgWithTheSameParamInd = symbols[token.str];
        							else {
        								symbols.insert(std::pair<std::string, int>(token.str, func.numberOfVectorArgs));
        								arg.firstVectorArgWithTheSameParamInd = func.numberOfVectorArgs;
        							}
        							++func.numberOfVectorArgs;
                                    break;
                                }
                                case E_NUMBER : {
                                    if (token.numType == E_REAL) throw 0; // !!!
                                    arg.varLength = false;
                                    std::stringstream ss(token.str);
                                    ss >> arg.length;
                                    ++func.numberOfVectorArgs;
                                    break;
                                }
                            }
							arg.type = E_ARGUMENT_TYPE_REAL_VECTOR;
						} else
							arg.type = E_ARGUMENT_TYPE_REAL;
						break;
					}
					case 'Z' : {
						if (GetNextToken(tokens, pos).str[0] == '^') {
							++pos;
							STokenInfo token = GetNextToken(tokens, pos);
							++pos;
							switch (token.type) {
                                case E_VARIABLE : {
                                    arg.varLength = true;
                                    if (symbols.find(token.str) != symbols.end()) arg.firstVectorArgWithTheSameParamInd = symbols[token.str];
        							else {
        								symbols.insert(std::pair<std::string, int>(token.str, func.numberOfVectorArgs));
        								arg.firstVectorArgWithTheSameParamInd = func.numberOfVectorArgs;
        							}
        							++func.numberOfVectorArgs;
                                    break;
                                }
                                case E_NUMBER : {
                                    if (token.numType != E_INTEGER) throw 0;
                                    arg.varLength = false;
                                    std::stringstream ss(token.str);
                                    ss >> arg.length;
                                    ++func.numberOfVectorArgs;
                                    break;
                                }
                            }
							arg.type = E_ARGUMENT_TYPE_INTEGER_VECTOR;
						} else
							arg.type = E_ARGUMENT_TYPE_INTEGER;
						break;
					}
					default : {
						throw 0;
					}
				}
    			token = GetNextToken(tokens, pos);
    			++pos;
    			if (token.token == E_TOWARDS) done = true;
    			else {
                    if (token.str != "x") throw 0;
                }
    			args.push_back(arg);
			} while (!done);
			
			token = GetNextToken(tokens, pos);
			++pos;
			switch (token.str[0]) {
				case 'R' : {
					if (GetNextToken(tokens, pos).str[0] == '^') {
						++pos;
						STokenInfo token = GetNextToken(tokens, pos);
						++pos;
						switch (token.type) {
                            case E_VARIABLE : {
                                func.retVal.varLength = true;
                                if (symbols.find(token.str) != symbols.end()) func.retVal.firstVectorArgWithTheSameParamInd = symbols[token.str];
    							else
    						        throw 0;
                                break;
                            }
                            case E_NUMBER : {
                                if (token.numType != E_INTEGER) throw 0;
                                func.retVal.varLength = false;
                                std::stringstream ss(token.str);
                                ss >> func.retVal.length;
                                break;
                            }
                        }
						func.retVal.type = E_ARGUMENT_TYPE_REAL_VECTOR;
					} else
						func.retVal.type = E_ARGUMENT_TYPE_REAL;
					break;
				}
				case 'Z' : {
					if (GetNextToken(tokens, pos).str[0] == '^') {
						++pos;
						STokenInfo token = GetNextToken(tokens, pos);
						++pos;
						switch (token.type) {
                            case E_VARIABLE : {
                                func.retVal.varLength = true;
                                if (symbols.find(token.str) != symbols.end()) func.retVal.firstVectorArgWithTheSameParamInd = symbols[token.str];
    							else
    						        throw 0;
                                break;
                            }
                            case E_NUMBER : {
                                if (token.numType != E_INTEGER) throw 0;
                                func.retVal.varLength = false;
                                std::stringstream ss(token.str);
                                ss >> func.retVal.length;
                                break;
                            }
                        }
						func.retVal.type = E_ARGUMENT_TYPE_INTEGER_VECTOR;
					} else
						func.retVal.type = E_ARGUMENT_TYPE_INTEGER;
					break;
				}
				default : {
					throw 0;
				}
            }
            
            // *** *** ***
            
            func.args = std::vector<SArgumentInfo>(args.begin(), args.end());
            
            // Creating the function signature.
            switch (func.retVal.type) {
                case E_ARGUMENT_TYPE_REAL : {
                    oss << "double ";
                    break;
                }
                case E_ARGUMENT_TYPE_REAL_VECTOR : {
                    oss << "vec ";
                    break;
                }
                case E_ARGUMENT_TYPE_INTEGER : {
                    oss << "int ";
                    break;
                }
                case E_ARGUMENT_TYPE_INTEGER_VECTOR : {
                    oss << "ivec ";
                    break;
                }
            }
            oss << func.name << "(";
            int varIndOld = aVarInd;
            for (int i = 0; i < func.args.size(); ++i) {
                int varNew = GetNextVar();
                if (i != 0) oss << ", ";
                switch (func.args[i].type) {
                    case E_ARGUMENT_TYPE_REAL : {
                        oss << "double $" << varNew;
                        break;
                    }
                    case E_ARGUMENT_TYPE_REAL_VECTOR : {
                        oss << "vec $" << varNew;
                        break;
                    }
                    case E_ARGUMENT_TYPE_INTEGER : {
                        oss << "int $" << varNew;
                        break;
                    }
                    case E_ARGUMENT_TYPE_INTEGER_VECTOR : {
                        oss << "ivec $" << varNew;
                        break;
                    }
                }
            }
            oss << ") {" << std::endl;
            
            // Determining the lengths of the vector arguments and checking if they match the domain. In the case of arguments with
            // the parameters we check if two arguments sharing the same parameter have the same length.
            int vectorArgNumber = 0;
            for (int i = 0; i < func.args.size(); ++i) {
                switch (func.args[i].type) {
                    case E_ARGUMENT_TYPE_REAL_VECTOR : {
                        int varNew = GetNextVar();
                        oss << std::string(tabSize, ' ') << "int $" << varNew << " = $" << varIndOld + i << ".n_elem;" << std::endl;
                        if (func.args[i].varLength) {
                            if (vectorArgNumber != func.args[i].firstVectorArgWithTheSameParamInd)
                                oss << std::string(tabSize, ' ') << "if ($" << varNew << " != $" << varIndOld + func.args.size() + func.args[i].firstVectorArgWithTheSameParamInd << ") throw 0;" << std::endl;
                        } else
                            oss << std::string(tabSize, ' ') << "if ($" << varNew << " != " << func.args[i].length << ") throw 0;" << std::endl;
                        ++vectorArgNumber;
                        break;
                    }
                    case E_ARGUMENT_TYPE_INTEGER_VECTOR : {
                        int varNew = GetNextVar();
                        oss << std::string(tabSize, ' ') << "int $" << varNew << " = $" << varIndOld + i << ".n_elem;" << std::endl;
                        if (func.args[i].varLength) {
                            if (vectorArgNumber != func.args[i].firstVectorArgWithTheSameParamInd)
                                oss << std::string(tabSize, ' ') << "if ($" << varNew << " != $" << varIndOld + func.args.size() + func.args[i].firstVectorArgWithTheSameParamInd << ") throw 0;" << std::endl;
                        } else
                            oss << std::string(tabSize, ' ') << "if ($" << varNew << " != " << func.args[i].length << ") throw 0;" << std::endl;
                        ++vectorArgNumber;
                        break;
                    }
                }
            }
            
            // Creating the variable for the return value.
            int retValInd = GetNextVar();
            switch (func.retVal.type) {
                case E_ARGUMENT_TYPE_REAL : {
                    oss << std::string(tabSize, ' ') << "double $" << retValInd << ";" << std::endl;
                    break;
                }
                case E_ARGUMENT_TYPE_REAL_VECTOR : {
                    if (func.retVal.varLength)
                        oss << std::string(tabSize, ' ') << "vec $" << retValInd << "($" << varIndOld + func.args.size() + func.retVal.firstVectorArgWithTheSameParamInd << ");" << std::endl;
                    else
                        oss << std::string(tabSize, ' ') << "vec $" << retValInd << "(" << func.retVal.length << ");" << std::endl;
                    break;
                }
                case E_ARGUMENT_TYPE_INTEGER : {
                    oss << std::string(tabSize, ' ') << "int $" << retValInd << ";" << std::endl;
                    break;
                }
                case E_ARGUMENT_TYPE_INTEGER_VECTOR : {
                    if (func.retVal.varLength)
                        oss << std::string(tabSize, ' ') << "ivec $" << retValInd << "($" << varIndOld + func.args.size() + func.retVal.firstVectorArgWithTheSameParamInd << ");" << std::endl;
                    else
                        oss << std::string(tabSize, ' ') << "ivec $" << retValInd << "(" << func.retVal.length << ");" << std::endl;
                    break;
                }
            }
            
            // *** *** ***
            
            // Switch with different versions of the function depending on the index value.
            if ((func.retVal.type == E_ARGUMENT_TYPE_REAL_VECTOR) || (func.retVal.type == E_ARGUMENT_TYPE_INTEGER_VECTOR)) {
                int varNew = GetNextVar();
                oss << std::string(tabSize, ' ') << "for (int $";
                oss << varNew << " = 1; $";
                if (func.retVal.varLength) oss << varNew << " <= $" << varIndOld + func.args.size() + func.retVal.firstVectorArgWithTheSameParamInd << "; ++$";
                else
                    oss << varNew << " <= " << func.retVal.length << "; ++$";
                oss << varNew << ") {" << std::endl;
                
                oss << std::string(tabSize * 2, ' ') << "switch ($" << varNew << ") {" << std::endl;
                
                bool done = false;
                bool indexVariableOccurred = false;
                do {
                    if (GetNextToken(tokens, pos).str != func.name) throw 0;
                    ++pos;
                    if (GetNextToken(tokens, pos).str[0] != '^') throw 0;
                    ++pos;
                    STokenInfo token = GetNextToken(tokens, pos);
                    ++pos;
                    switch (token.type) {
                        case E_NUMBER : {
                            // The function definition with the index variable must be the last definition.
                            if (indexVariableOccurred) throw 0;
                            if (token.numType == E_REAL) throw 0;
                            std::stringstream ss(token.str);
                            int val;
                            ss >> val;
                            oss << std::string(tabSize * 3, ' ') << "case " << val << " : {" << std::endl;
                            
                            // !!! !!! !!!
                            // !!! !!! !!!
                            bool done = false;
                            do {
                                if (GetNextToken(tokens, pos).token == E_RIGHT_PAR) done = true;
                                ++pos;
                            } while (!done);
                            // !!! !!! !!!
                            // !!! !!! !!!
                            
                            if (GetNextToken(tokens, pos).token != E_ASSIGN) throw 0;
                            ++pos;
                            
                            std::string code;
                            int varRes;
                            if (!ParseExpression(tokens, pos, code, varRes)) throw 0;
                            oss << MakeIndent(code, ' ', tabSize * 4);
                            
                            oss << std::string(tabSize * 4, ' ') << "$" << retValInd << "(" << "$" << varNew << " - 1) = $" << varRes << ";" << std::endl;
                            oss << std::string(tabSize * 4, ' ') << "break;" << std::endl;
                            oss << std::string(tabSize * 3, ' ') << "}" << std::endl;
                            break;
                        }
                        case E_VARIABLE : {
                            // The function definition with the index variable cannot occur more than one time.
                            if (indexVariableOccurred) throw 0;
                            oss << std::string(tabSize * 3, ' ') << "default : {" << std::endl;
                            
                            // !!! !!! !!!
                            // !!! !!! !!!
                            // !!! !!! !!!
                            bool done = false;
                            do {
                                if (GetNextToken(tokens, pos).token == E_RIGHT_PAR) done = true;
                                ++pos;
                            } while (!done);
                            // !!! !!! !!!
                            // !!! !!! !!!
                            // !!! !!! !!!
                            
                            if (GetNextToken(tokens, pos).token != E_ASSIGN) throw 0;
                            ++pos;
                            
                            std::string code;
                            int varRes;
                            if (!ParseExpression(tokens, pos, code, varRes)) throw 0;
                            oss << MakeIndent(code, ' ', tabSize * 4);
                            
                            oss << std::string(tabSize * 4, ' ') << "$" << retValInd << "(" << "$" << varNew << " - 1) = $" << varRes << ";" << std::endl;
                            oss << std::string(tabSize * 3, ' ') << "}" << std::endl;
                            indexVariableOccurred = true;
                            break;
                        }
                        default : {
                            throw 0;
                        }
                    }
                    
                    if (pos < tokens.size()) {
                        if (GetNextToken(tokens, pos).str != func.name) done = true;
                    } else
                        done = true;
                } while (!done);
                
                oss << std::string(tabSize * 2, ' ') << "}" << std::endl;
                oss << std::string(tabSize, ' ') << "}" << std::endl;    
            }
            
            // *** *** ***
            
            // Creating the function ending.
            oss << std::string(tabSize, ' ') << "return $" << retValInd << ";" << std::endl;
            oss << "}" << std::endl;
            oss << std::endl;
            
            //break; // !!!
        }
        code = oss.str();      
    } catch (...) {
        return false;
    }
    return true;
}