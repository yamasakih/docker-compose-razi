from rdkit import Chem

from razi.rdkit_postgresql.types import Mol, Bfp
from razi.rdkit_postgresql.functions import atompairbv_fp, torsionbv_fp, morganbv_fp, featmorganbv_fp

from sqlalchemy import create_engine, Column, Index, Integer, String
from sqlalchemy.ext.declarative import declarative_base


engine = create_engine('postgresql://postgres@db:5432/souyakuRB2018')
Base = declarative_base(bind=engine)


class Compound(Base):
    __tablename__ = 'compounds'
    
    id = Column(Integer, primary_key=True)
    structure = Column(Mol)
    name = Column(String)
    library_id = Column(Integer)
    mfp2 = Column(Bfp)
    ffp2 = Column(Bfp)
    atompair = Column(Bfp)
    torsion = Column(Bfp)

    __table_args__ = (
        Index('structure_idx', 'structure',
                    postgresql_using='gist'),
    )
    
    def __init__(self, structure, name=None, library_id=None):
        self.name = name
        self.library_id = library_id
        if isinstance(structure, Chem.Mol):
            self.structure = Chem.MolToSmiles(structure)
        elif isinstance(structure, str):
            self.structure = structure
        self.mfp2 = morganbv_fp(self.structure, 2)
        self.ffp2 = featmorganbv_fp(self.structure, 2)
        self.atompair = atompairbv_fp(self.structure)
        self.torsion = torsionbv_fp(self.structure)
        
    def __repr__(self):
        if isinstance(self.structure, Chem.Mol):
            return '(%s) < %s >' % (self.name, Chem.MolToSmiles(self.structure))
        return '(%s) < %s >' % (self.name, self.structure)


class Library(Base):
    __tablename__ = 'libraries'
    
    id = Column(Integer, primary_key=True)
    name = Column(String)
    
    def __init__(self, id, name):
        self.id = id
        self.name = name
        
    def __repr__(self):
        return f'{self.name}'


if __name__ == '__main__':
    Base.metadata.create_all()
