version 1.0

struct Person {
  String name
  Int age
}

task printPerson {
  input {
    Person person
    Int years
  }

  command {
    echo "hello my name is ${person.name} and I am ${person.age} years old, although writing WDL has aged me ~{years} years"
  }

  output {
    Person older_person = object {
      name: person.name,
      age: person.age + years
    }
    String message = read_string(stdout())
  }
}

workflow wf_person {
  input {
    Array[Person] people
    Array[Int] years
  }

  scatter (item in zip(people, years)) {
    call printPerson {
      input:
        person = item.left,
        years =  item.right
    }
  }

  output {
    Array[Person] older_people = printPerson.older_person
    Array[String] messages = printPerson.message
  }
}
