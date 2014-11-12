// $Id$
#ifndef INCLUDED_CONSERVATIVE_MAPPING_H__
#define INCLUDED_CONSERVATIVE_MAPPING_H__

#include "AbuttingInterface.h"

class ConservativeMapping
{
public:
    ConservativeMapping(
        const AbuttingInterface& interface
        );

protected:

private:
    const AbuttingInterface& mInterface;
};

#endif // INCLUDED_CONSERVATIVE_MAPPING_H__

