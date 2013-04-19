forAll(PhiBoundaryTypes, patchi)
{
    if (PhiBoundaryTypes[patchi] == "fixedGradient")
    {
        fixedGradientFvPatchField<scalar>& pf =
            const_cast<fixedGradientFvPatchField<scalar>&>
            (
                refCast<const fixedGradientFvPatchField<scalar> >
                (
                    Phi.boundaryField()[patchi]
                )
            );
        pf.gradient() = mesh.boundary()[patchi].nf() & U.boundaryField()[patchi];
    }
}
