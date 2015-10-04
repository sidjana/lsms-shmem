#include <string.h>
#include "Real.hpp"
#include "LuaSupport.hpp"

bool luaGetStrN(lua_State *L, const char *name, char *s, int n)
{
  lua_getglobal(L,name);
  if(!lua_isstring(L,-1)) {lua_pop(L,1); return false;}
  strncpy(s, lua_tostring(L,-1), n);
  lua_pop(L,1);
  return true;
}

bool luaGetReal(lua_State *L, const char *name, Real *val)
{
  lua_getglobal(L,name);
  if(!lua_isnumber(L,-1)) {lua_pop(L,1); return false;}
  *val=lua_tonumber(L,-1);
  lua_pop(L,1);
  return true;
}

bool luaGetInteger(lua_State *L, const char *name, int *val)
{
  lua_getglobal(L,name);
  if(!lua_isnumber(L,-1)) {lua_pop(L,1); return false;}
  *val=lua_tointeger(L,-1);
  lua_pop(L,1);
  return true;
}

bool luaGetBoolean(lua_State *L, const char *name, bool *val)
{
  lua_getglobal(L,name);
  if(!lua_isboolean(L,-1)) {lua_pop(L,1); return false;}
  *val=lua_toboolean(L,-1);
  lua_pop(L,1);
  return true;
}

bool luaGetFieldInTable(lua_State *L, const char *name, const char *field)
{
  lua_getglobal(L,name);
  lua_pushstring(L,field);
  lua_gettable(L,-2);
  return true;
}

bool luaGetFieldFromStack(lua_State *L, const char *field)
{
  lua_pushstring(L,field);
  lua_gettable(L,-2);
  return true;
}

bool luaGetPositionInTable(lua_State *L, const char *name, int idx)
{
  lua_getglobal(L,name);
  lua_pushinteger(L,idx);
  lua_gettable(L,-2);
  return true;
}

bool luaGetRealFieldInTable(lua_State *L, const char *name, const char *field, Real *val)
{
  luaGetFieldInTable(L,name,field);
  if(!lua_isnumber(L,-1)) {lua_pop(L,1); return false;}
  *val=lua_tonumber(L,-1);
  lua_pop(L,1);
  return true;
}

bool luaGetIntegerFieldInTable(lua_State *L, const char *name, const char *field, int *val)
{
  luaGetFieldInTable(L,name,field);
  if(!lua_isnumber(L,-1)) {lua_pop(L,1); return false;}
  *val=lua_tointeger(L,-1);
  lua_pop(L,1);
  return true;
}

bool luaGetRealPositionInTable(lua_State *L, const char *name, int idx, Real *val)
{
  luaGetPositionInTable(L,name,idx);
  if(!lua_isnumber(L,-1)) {lua_pop(L,1); return false;}
  *val=lua_tonumber(L,-1);
  lua_pop(L,1);
  return true;
}

bool luaGetRealPositionFromStack(lua_State *L, int idx, Real *val)
{
  lua_pushinteger(L,idx);
  lua_gettable(L,-2);
  if(!lua_isnumber(L,-1)) {lua_pop(L,1); return false;}
  *val=lua_tonumber(L,-1);
  lua_pop(L,1);
  return true;
}

bool luaGetStrNFromStack(lua_State *L, const char *name, char *s, int n)
{
  lua_pushstring(L,name);
  lua_gettable(L,-2);
  if(!lua_isstring(L,-1)) {lua_pop(L,1); return false;}
  strncpy(s, lua_tostring(L,-1), n);
  lua_pop(L,1);
  return true;
}

bool luaGetRealFieldFromStack(lua_State *L, const char *field, Real *val)
{
  lua_pushstring(L,field);
  lua_gettable(L,-2);
  if(!lua_isnumber(L,-1)) {lua_pop(L,1); return false;}
  *val=lua_tonumber(L,-1);
  lua_pop(L,1);
  return true;
}

bool luaGetIntegerFieldFromStack(lua_State *L, const char *field, int *val)
{
  lua_pushstring(L,field);
  lua_gettable(L,-2);
  if(!lua_isnumber(L,-1)) {lua_pop(L,1); return false;}
  *val=lua_tointeger(L,-1);
  lua_pop(L,1);
  return true;
}
