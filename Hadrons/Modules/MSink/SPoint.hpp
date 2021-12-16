/*
 * SPoint.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Tristan Ueding <trisueding@web.de>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution
 * directory.
 */

/*  END LEGAL */

#ifndef Hadrons_MSink_SPoint_hpp_
#define Hadrons_MSink_SPoint_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

// The only difference to Module Point of MSink is that here the slicing in the end
// is done in z direction and not in the temporal direction and Point --> ZPoint

/******************************************************************************
 *                                   ZPoint                                    *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSink)

class ZPointPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ZPointPar,
                                    std::string, mom);
};

template <typename FImpl>
class ZPoint: public Module<ZPointPar>
{
public:
    BASIC_TYPE_ALIASES(FImpl,);
    SINK_TYPE_ALIASES();
public:
    // constructor
    ZPoint(const std::string name);
    // destructor
    virtual ~ZPoint(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    bool        hasPhase_{false};
    std::string momphName_;
};

MODULE_REGISTER_TMP(SPoint,       ZPoint<FIMPL>,        MSink);
MODULE_REGISTER_TMP(ScalarSPoint, ZPoint<ScalarImplCR>, MSink);

/******************************************************************************
 *                          ZPoint implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
ZPoint<FImpl>::ZPoint(const std::string name)
: Module<ZPointPar>(name)
, momphName_ (name + "_momph")
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> ZPoint<FImpl>::getInput(void)
{
    std::vector<std::string> in;

    return in;
}

template <typename FImpl>
std::vector<std::string> ZPoint<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void ZPoint<FImpl>::setup(void)
{
    envTmpLat(LatticeComplex, "coor");
    envCacheLat(LatticeComplex, momphName_);
    envCreate(SinkFn, getName(), 1, nullptr);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void ZPoint<FImpl>::execute(void)
{
    LOG(Message) << "Setting up zpoint sink function for momentum ["
                 << par().mom << "]" << std::endl;

    auto &ph = envGet(LatticeComplex, momphName_);

    if (!hasPhase_)
    {
        Complex           i(0.0,1.0);
        std::vector<Real> p;

        envGetTmp(LatticeComplex, coor);
        p  = strToVec<Real>(par().mom);
        ph = Zero();
        for(unsigned int mu = 0; mu < p.size(); mu++)
        {
            // In the last iteration take the temporal component
            // since this source is for a spatial correlator
            if(mu == 2){
                mu+=1;
            }
            LatticeCoordinate(coor, mu);
            ph = ph + (p[mu]/env().getDim(mu))*coor;
        }
        ph = exp((Real)(2*M_PI)*i*ph);
        hasPhase_ = true;
    }
    auto sink = [&ph](const PropagatorField &field)
    {
        SlicedPropagator res;
        PropagatorField  tmp = ph*field;

        sliceSum(tmp, res, Zp);

        return res;
    };
    envGet(SinkFn, getName()) = sink;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSink_SPoint_hpp_
