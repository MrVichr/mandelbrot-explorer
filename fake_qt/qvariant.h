#pragma once
#include "QObject"
#include "QString"
#include "QVector"
#include <variant>
#include <cstdlib>

// Minimal QVariant stub — stores int, double, QString, or bool
class QVariant {
    enum class Type { Null, Int, Double, String, Bool };
    Type _type = Type::Null;
    std::variant<std::monostate, int, double, QString, bool> _val;

public:
    QVariant() : _type(Type::Null), _val(std::monostate{}) {}
    QVariant(int v)           : _type(Type::Int),    _val(v) {}
    QVariant(double v)        : _type(Type::Double),  _val(v) {}
    QVariant(float v)         : _type(Type::Double),  _val((double)v) {}
    QVariant(const QString &v): _type(Type::String),  _val(v) {}
    QVariant(const char *v)   : _type(Type::String),  _val(QString(v)) {}
    QVariant(bool v)          : _type(Type::Bool),    _val(v) {}

    bool isNull()  const { return _type == Type::Null; }
    bool isValid() const { return _type != Type::Null; }

    double toDouble(bool *ok = nullptr) const {
        if (ok) *ok = true;
        switch (_type) {
            case Type::Double: return std::get<double>(_val);
            case Type::Int:    return (double)std::get<int>(_val);
            case Type::Bool:   return std::get<bool>(_val) ? 1.0 : 0.0;
            case Type::String: {
                const std::string &s = std::get<QString>(_val).toStdString();
                try { return std::stod(s); }
                catch (...) { if (ok) *ok = false; return 0.0; }
            }
            default: if (ok) *ok = false; return 0.0;
        }
    }

    int toInt(bool *ok = nullptr) const {
        if (ok) *ok = true;
        switch (_type) {
            case Type::Int:    return std::get<int>(_val);
            case Type::Double: return (int)std::get<double>(_val);
            case Type::Bool:   return std::get<bool>(_val) ? 1 : 0;
            case Type::String: {
                const std::string &s = std::get<QString>(_val).toStdString();
                try { return std::stoi(s); }
                catch (...) { if (ok) *ok = false; return 0; }
            }
            default: if (ok) *ok = false; return 0;
        }
    }

    QString toString() const {
        switch (_type) {
            case Type::String: return std::get<QString>(_val);
            case Type::Int:    return QString::number(std::get<int>(_val));
            case Type::Double: return QString::number(std::get<double>(_val));
            case Type::Bool:   return std::get<bool>(_val) ? QString("true") : QString("false");
            default: return {};
        }
    }

    bool toBool() const {
        switch (_type) {
            case Type::Bool:   return std::get<bool>(_val);
            case Type::Int:    return std::get<int>(_val) != 0;
            case Type::Double: return std::get<double>(_val) != 0.0;
            case Type::String: {
                const QString &s = std::get<QString>(_val);
                return s == "true" || s == "1";
            }
            default: return false;
        }
    }
};

// QVariantMap = QMap<QString, QVariant>
using QVariantMap = QMap<QString, QVariant>;
