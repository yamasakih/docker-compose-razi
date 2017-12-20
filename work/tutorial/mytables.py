from rdkit import Chem

from razi.rdkit_postgresql.types import Mol, Bfp
from razi.rdkit_postgresql.functions import atompairbv_fp, torsionbv_fp, morganbv_fp

from sqlalchemy import create_engine, Column, Index, Integer, String
from sqlalchemy.ext.declarative import declarative_base


engine = create_engine('postgresql://postgres@db:5432/postgres')
Base = declarative_base(bind=engine)


class Compound(Base):
    __tablename__ = 'compounds'
    
    id = Column(Integer, primary_key=True)
    name = Column(String)
    structure = Column(Mol)
    atompair = Column(Bfp)
    torsion = Column(Bfp)
    morgan = Column(Bfp)
    
    __table_args__ = (
        Index('compounds_structure', 'structure',
                    postgresql_using='gist'),
        )
    
    def __init__(self, name, structure):
        self.name = name
        if isinstance(structure, Chem.Mol):
            self.structure = Chem.MolToSmiles(structure)
        elif isinstance(structure, str):
            self.structure = structure
        self.atompair = atompairbv_fp(self.structure)
        self.torsion = torsionbv_fp(self.structure)
        self.morgan = morganbv_fp(self.structure, 2)
        
    def __repr__(self):
        if isinstance(self.structure, Chem.Mol):
            return '(%s) < %s >' % (self.name, Chem.MolToSmiles(self.structure))
        return '(%s) < %s >' % (self.name, self.structure)


if __name__ == '__main__':
    Base.metadata.create_all()
