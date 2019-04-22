#include "Functions.h"

Functions::Functions()
{
    //ctor
}

// char Global::charFromFunction(int f)
//{
//  switch(f)
//    {
//    case Global::Plus:
//      return '+';
//    case Global::Minus:
//      return '-';
//    case Global::Times:
//      return '*';
//    case Global::Divide:
//      return '/';
//    case Global::Sqrt:
//      return 'Q';
//    case Global::Exp:
//      return 'E';
//    case Global::Sin:
//      return 'S';
//    case Global::LessThan:
//      return 'L';
//    case Global::GreaterThan:
//      return 'G';
//    case Global::Abs:
//      return 'A';
//    case Global::Log:
//      return 'O';
//    default:
//      return '!';
//    }
//}
//
//char Global::charFromInt(int i)
//{
//  char c = '!';
//  if(i < Constant ) {
//    c = Global::charFromFunction(i);
//  }
//  else if( i == Constant ) {
//    c = '?';
//  }
//  else if( i >= Variable ) {
//    int ic = (i-Variable) % 26; // we map them all to a-z
//    c = 'a' + ic;
//  }
//  return c;
//}
