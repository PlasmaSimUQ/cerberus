#ifndef MFP_FACTORY_H
#define MFP_FACTORY_H

#include "sol.hpp"

#include <functional>
#include <map>
#include <memory>
#include <string>

template <class T> class ClassFactory
{
  public:
    typedef std::function<std::unique_ptr<T>(const sol::table&)> Builder;

    /// returns true if the registration succeeded, false otherwise
    bool Register(std::string const& key, Builder const& builder)
    {
        _map.insert(std::make_pair(key, builder));
        return true;
    }

    /// returns a pointer to a new instance
    std::unique_ptr<T> Build(const std::string& key, const sol::table& def = sol::table()) const
    {
        auto range = _map.equal_range(key);
        for (auto it = range.first; it != range.second; ++it) {
            std::unique_ptr<T> obj = it->second(def);
            if (obj) return obj;
        }

        return 0;
    }

    const std::multimap<std::string, Builder>& getRegistered() const { return _map; }

    const std::vector<std::string> getKeys() const
    {
        std::vector<std::string> keys;
        for (const auto& m : _map) { keys.push_back(m.first); }
        return keys;
    }

  private:
    std::multimap<std::string, Builder> _map;
};

#endif  // MFP_FACTORY_H
