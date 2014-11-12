// $Id: PatchExchanger.h 39 2010-06-21 02:42:48Z kato $
#ifndef INCLUDED_PATCH_EXCHANGER_H__
#define INCLUDED_PATCH_EXCHANGER_H__

class PatchExchanger
{
public:
    virtual ~PatchExchanger() {}

    virtual void Start() = 0;
    virtual void Finish() = 0;

protected:

private:
};

#endif // INCLUDED_PATCH_EXCHANGER_H__

