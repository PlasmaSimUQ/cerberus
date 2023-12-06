#include "LuaContext.H"

bool lua_variable_exists(LuaContext& lua, const std::vector<std::string>& path, std::string& loc)
{
    if (path.empty()) { return false; }

    bool exists = false;
    loc = "";
    for (int i = 0; i < path.size(); ++i) {
        if (i == 0) {
            loc += path[i];
        } else {
            loc += "['" + path[i] + "']";
        }

        exists = lua.executeCode<int>("if not " + loc + " then return 0 else return 1 end");
        if (!exists) break;
    }

    return exists;
}
