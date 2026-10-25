// moc_signals.cpp — signal method bodies that Qt's moc would normally generate.
// Each signal body simply calls QObject::emitSignal(), which looks up all
// connections registered for (this, PMF) and dispatches them.

#include "../MandelEvaluator.hpp"
#include "../MandelModel.hpp"
#include "../JuliaModel.hpp"
#include "../LaguerreModel.hpp"

// ── MandelEvaluatorThread signals ──────────────────────────────────────────

void MandelEvaluatorThread::doneMandelThreaded(MandelEvaluatorThread *me)
{ QObject::emitSignal(this, &MandelEvaluatorThread::doneMandelThreaded, me); }

void MandelEvaluatorThread::doneJuliaThreaded(MandelEvaluatorThread *me)
{ QObject::emitSignal(this, &MandelEvaluatorThread::doneJuliaThreaded, me); }

void MandelEvaluatorThread::doneLaguerreThreaded(MandelEvaluatorThread *me)
{ QObject::emitSignal(this, &MandelEvaluatorThread::doneLaguerreThreaded, me); }

// ── MandelModel signals ────────────────────────────────────────────────────

void MandelModel::selectedPaintStyleChanged()
{ QObject::emitSignal(this, &MandelModel::selectedPaintStyleChanged); }

void MandelModel::selectedPrecisionChange()
{ QObject::emitSignal(this, &MandelModel::selectedPrecisionChange); }

void MandelModel::extAngleZoomChange()
{ QObject::emitSignal(this, &MandelModel::extAngleZoomChange); }

void MandelModel::triggerComputeThreaded(int epoch)
{ QObject::emitSignal(this, &MandelModel::triggerComputeThreaded, epoch); }

// ── JuliaModel signals ─────────────────────────────────────────────────────

void JuliaModel::selectedPaintStyleChanged()
{ QObject::emitSignal(this, &JuliaModel::selectedPaintStyleChanged); }

void JuliaModel::selectedPrecisionChange()
{ QObject::emitSignal(this, &JuliaModel::selectedPrecisionChange); }

void JuliaModel::selectedEdePatchAlgoChange()
{ QObject::emitSignal(this, &JuliaModel::selectedEdePatchAlgoChange); }

void JuliaModel::selectedExteriorColoringChange()
{ QObject::emitSignal(this, &JuliaModel::selectedExteriorColoringChange); }

void JuliaModel::orbit_frozenChange()
{ QObject::emitSignal(this, &JuliaModel::orbit_frozenChange); }

void JuliaModel::extAngleZoomChange()
{ QObject::emitSignal(this, &JuliaModel::extAngleZoomChange); }

void JuliaModel::triggerJuliaThreaded(int epoch, int juliaPeriod)
{ QObject::emitSignal(this, &JuliaModel::triggerJuliaThreaded, epoch, juliaPeriod); }

// ── LaguerreModel signals ──────────────────────────────────────────────────

void LaguerreModel::selectedPaintStyleChanged()
{ QObject::emitSignal(this, &LaguerreModel::selectedPaintStyleChanged); }

void LaguerreModel::selectedPrecisionChange()
{ QObject::emitSignal(this, &LaguerreModel::selectedPrecisionChange); }

void LaguerreModel::triggerLaguerreThreaded(int epoch)
{ QObject::emitSignal(this, &LaguerreModel::triggerLaguerreThreaded, epoch); }
